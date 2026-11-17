"""Tiempos por fase, impresión de resultados y fila del CSV común."""

import os
from dataclasses import dataclass, fields

COLUMNAS_CSV = ("version,distribucion,N,pasos,theta,epsilon,dt,trabajadores,"
                "t_construccion,t_ordenamiento,t_fuerzas,t_integracion,t_total,"
                "energia_inicial,energia_final,error_energia")


@dataclass
class TiemposPorFase:
    """Segundos por fase. En un paso son los de ese paso; en el resumen, la suma."""
    construccion: float = 0.0
    ordenamiento: float = 0.0
    fuerzas: float = 0.0
    integracion: float = 0.0
    total: float = 0.0
    desbalance: float = 0.0           # trabajador más lento / promedio (1 = perfecto)
    construccion_local: float = 0.0   # estrategia local: incluida en fuerzas

    def sumar(self, otro):
        for campo in fields(self):
            setattr(self, campo.name, getattr(self, campo.name) + getattr(otro, campo.name))


def error_relativo(valor, referencia):
    return abs(valor - referencia) / abs(referencia)


def imprimir_entorno(version, config, num_particulas):
    print(f"Barnes-Hut, versión {version}")
    print(f"  entrada:             {config.entrada} (distribución {config.distribucion})")
    print(f"  N:                   {num_particulas}")
    print(f"  pasos:               {config.pasos}")
    print(f"  θ:                   {config.theta:g}")
    print(f"  ε:                   {config.epsilon:g}")
    print(f"  Δt:                  {config.dt:g}")
    print(f"  trabajadores:        {config.trabajadores}")
    print(f"  estrategia del árbol: {config.estrategia}")
    print(f"  reparto:             {config.reparto}")
    print(f"  tamaño de grupo:     {config.tamano_grupo}")
    print(f"  núcleos disponibles: {os.cpu_count()}")
    print()


def imprimir_encabezado_pasos():
    print("Tiempos por paso (ms); desbalance = trabajador más lento / promedio en fuerzas")
    print(f"{'paso':>6} {'construcción':>13} {'ordenamiento':>13} {'fuerzas':>10} "
          f"{'integración':>12} {'total':>10} {'desbalance':>11}")


def imprimir_paso(paso, tiempos):
    print(f"{paso:>6} {1e3 * tiempos.construccion:>13.2f} {1e3 * tiempos.ordenamiento:>13.2f} "
          f"{1e3 * tiempos.fuerzas:>10.2f} {1e3 * tiempos.integracion:>12.2f} "
          f"{1e3 * tiempos.total:>10.2f} {tiempos.desbalance:>11.2f}")


def imprimir_energia(paso, energia, energia_inicial):
    print(f"       energía en el paso {paso}: {energia:.10f} "
          f"(error relativo {error_relativo(energia, energia_inicial):.3e})")


def imprimir_resumen(acumulado, pasos, energia_inicial, energia_final, estrategia):
    print(f"\nTiempos acumulados de {pasos} pasos")
    print(f"  {'fase':<24} {'total (s)':>10} {'por paso (ms)':>14} {'porcentaje':>13}")
    for nombre, segundos in (("construcción del árbol", acumulado.construccion),
                             ("ordenamiento", acumulado.ordenamiento),
                             ("fuerzas", acumulado.fuerzas),
                             ("integración", acumulado.integracion),
                             ("total", acumulado.total)):
        por_paso = 1e3 * segundos / pasos if pasos > 0 else 0.0
        porcentaje = 100 * segundos / acumulado.total if acumulado.total > 0 else 0.0
        print(f"  {nombre:<24} {segundos:>10.3f} {por_paso:>14.2f} {porcentaje:>11.1f} %")
    if pasos > 0:
        print(f"  desbalance promedio en fuerzas: {acumulado.desbalance / pasos:.2f}")
        if estrategia == "local":
            print(f"  construcción del árbol en cada trabajador (incluida en fuerzas): "
                  f"{1e3 * acumulado.construccion_local / pasos:.2f} ms por paso")

    print("\nEnergía (calculada con el árbol)")
    print(f"  inicial:        {energia_inicial:.10f}")
    print(f"  final:          {energia_final:.10f}")
    print(f"  error relativo: {error_relativo(energia_final, energia_inicial):.3e}")


def agregar_fila_csv(config, version, num_particulas, acumulado, energia_inicial, energia_final):
    try:
        with open(config.csv, "a", encoding="utf-8") as archivo:
            # El encabezado solo va si el archivo está vacío, así varias
            # corridas de distintas versiones quedan en la misma tabla.
            if archivo.tell() == 0:
                archivo.write(COLUMNAS_CSV + "\n")
            archivo.write(
                f"{version},{config.distribucion},{num_particulas},{config.pasos},"
                f"{config.theta:g},{config.epsilon:g},{config.dt:g},{config.trabajadores},"
                f"{acumulado.construccion:.6f},{acumulado.ordenamiento:.6f},{acumulado.fuerzas:.6f},"
                f"{acumulado.integracion:.6f},{acumulado.total:.6f},"
                f"{energia_inicial:.10g},{energia_final:.10g},"
                f"{error_relativo(energia_final, energia_inicial):.6e}\n")
    except OSError:
        print(f"Aviso: no se pudo abrir {config.csv}; la medición no se guardó")
        return
    print(f"\nMedición agregada a {config.csv}")
