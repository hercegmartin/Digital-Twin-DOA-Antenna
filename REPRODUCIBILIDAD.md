# Arquitectura Computacional y Reproducibilidad

Este documento detalla los requerimientos de software y la estructura del flujo de datos necesarios para ejecutar el Gemelo Digital de forma local. El código fue diseñado priorizando la estabilidad del motor electromagnético y la automatización del post-proceso.

## 1. Entorno y Dependencias

El simulador FDTD no es una librería estándar de Python; requiere una instalación específica del motor **openEMS** en el sistema operativo.

### 1.1 Motor openEMS (Core Electromagnético)
* Es obligatorio tener instalado openEMS (binarios de Windows/Linux).
* El script inyecta la ruta de openEMS directamente en el `PATH` de Windows durante la ejecución para evitar conflictos de variables de entorno. 
* **Ruta por defecto en el script:** `C:\Users\herce\openEMS\openEMS`. *(Modificar la variable `openems_path` en la Sección 1 si tu instalación está en otro directorio).*

### 1.2 Entorno Python (Conda)
Se recomienda aislar las dependencias utilizando un entorno de Conda. Las librerías requeridas para el post-proceso DSP y la generación de métricas son estándar:
* `numpy`: Para el álgebra matricial, cálculo de SWR y algoritmos DOA (Descomposición en autovalores).
* `matplotlib`: Para la exportación silenciosa de los diagramas de radiación y el espectro MUSIC/Bartlett.
* `CSXCAD`: Interfaz de geometría (incluida con openEMS).

## 2. Flujo de Ejecución (Pipeline)

El script opera como una "caja negra" evaluadora. No requiere intervención manual una vez iniciado.

1. **Pre-procesamiento (Geometría):** Genera la cruz planar basándose en el espacio de parámetros provisto (frecuencia, dimensiones, materiales).
2. **Excitación Secuencial:** Lanza 5 simulaciones FDTD aisladas. En cada una, inyecta un pulso gaussiano en un puerto mientras los otros cuatro actúan como cargas adaptadas pasivas de $50 \Omega$.
3. **Validación Multi-fidelidad:** Antes de simular todo el barrido, el flag `verificar_primera_geometria = True` abre automáticamente la interfaz 3D (AppCSXCAD) en la primera iteración. Esto permite auditar visualmente el *mesh snapping* y la estructura física antes de comprometer horas de cómputo.

![Malla 3D y Geometría en AppCSXCAD](assets/mesh_3d.png)

4. **Transformación NF2FF:** Se proyectan los campos calculados a la zona de Fraunhofer, extrayendo el *Manifold* complejo.
5. **Post-proceso DSP:** Se inyecta ruido blanco (SNR definido por el usuario) y se calculan los espectros de Bartlett y MUSIC espacial.

## 3. Estructura de Salida (Outputs)

Para evitar la saturación de memoria en barridos largos, la generación de gráficos (`plt.show()`) está suprimida. Todos los resultados se escriben de forma no volátil en el disco.

Por cada simulación, el script genera una carpeta aislada que contiene:
* `simp_patch_array.xml`: El archivo de geometría crudo para inspección en AppCSXCAD.
* `1_S11_Central.png`: Gráfico de adaptación de impedancia y resonancia.
* `2_FoV_Fisico.png`: Corte del lóbulo de radiación (Planos E y H).
* `3_Espectro_DOA.png`: Superposición del pseudo-espectro de los algoritmos DOA.

Finalmente, el script unifica las métricas de todas las iteraciones en un archivo maestro `Reporte_Barrido.txt` con formato de tabla estándar, ideal para ser parseado posteriormente por un script de Optimización Bayesiana.