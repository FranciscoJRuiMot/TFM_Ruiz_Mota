# Promoting stability and resilience in mutualistic communities of evolving digital organisms

Código del Trabajo Fin de Máster (Máster en Análisis de Datos Ómicos y Biología de Sistemas, Universidad de Sevilla y Universidad Internacional de Andalucía).

## Contexto

En las comunidades microbianas, muchas especies viven de los metabolitos que liberan otras (*cross feeding*). En este trabajo estudiamos qué hace estables estas redes de intercambio y cómo cambian con la evolución, usando organismos digitales: programas que se replican, mutan y evolucionan en la plataforma [Avida](https://github.com/devosoft/avida). Cada especie realiza una combinación de dos funciones booleanas, que equivalen a reacciones metabólicas, y libera subproductos que consumen otras especies.

El trabajo tiene dos partes:

- **Competencia ecológica.** 1.000 redes de 10 especies simuladas sin mutación, con y sin competencia por recursos (fusionando el recurso de dos especies), para identificar qué propiedades de la red y de los organismos determinan la estabilidad de la comunidad.
- **Evolución.** 10 comunidades ancestrales de 9 especies, extraídas de un universo de 72 especies posibles, que evolucionan durante 50.000 *updates* en 10 réplicas cada una, para ver qué especies nuevas surgen y cómo cambia la red.

## Estructura del repositorio

```
├── Competencia_ecologica/   Experimentos de la parte ecológica
├── Evolucion/               Experimentos de la parte evolutiva
├── Analisis_estadistico/    Análisis de los resultados de las simulaciones
└── Figuras/                 Figuras del manuscrito
```

**`Competencia_ecologica/`** contiene el código que genera las redes de interacción y ejecuta en Avida las simulaciones con y sin competencia.

**`Evolucion/`** contiene el código que genera el universo de interacciones, construye las comunidades ancestrales, ejecuta en Avida las réplicas evolutivas y analiza los fenotipos resultantes.

**`Analisis_estadistico/`** contiene un script por cada parte del trabajo. Cada uno reproduce las cifras y tablas del apartado de resultados correspondiente:

| Script | Manuscrito | Salidas |
|---|---|---|
| `analisis_ecologia.R` | Apartado 3.1 | Tabla 1, Tabla S2 |
| `analisis_evolucion.R` | Apartado 3.2 | Tablas 2 y 3, Tabla S1 |

Las cifras se muestran en consola, ordenadas por subsección, y se guardan en `resultados/` junto con las tablas en formato CSV. Los datos de entrada de ambos scripts están en `datos/`.

**`Figuras/`** contiene `figuras.R`, que genera las Figuras 1–5 y S1–S2.