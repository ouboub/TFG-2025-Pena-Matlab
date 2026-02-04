# Cómo usar el repositorio

Este repositorio contiene implementaciones de distintos métodos numéricos para resolver la ecuación de Lane-Emden. No existe un programa principal, cada método se ejecuta de forma separada, salvo en el caso de los métodos SLM y QLM, que se prueban conjuntamente.


## Métodos SLM y QLM

Los métodos SLM y QLM se prueban mediante el script:

```
> laneemdenprueba
```

Este archivo `laneemdenprueba.m` está diseñado para:

-   Ejecutar el método **SLM**.

-   Ejecutar el método **QLM**.

No incluye los métodos `ode45` ni `RK4`. Su función es únicamente validar y analizar el comportamiento de los métodos iterativos. Cuando se ejecuta este script se obtiene como resultado una tabla con la solución aproximada para los distintos valores del índice politrópico en cada iteración. Depende del método que se quiera probar, se deberá de introducir `laneemdenSLM.m` o `laneemdenNEWTONKANTOROVICH.m` en la línea correcta del script.


## Métodos clásicos `ode45` y `RK4`

Los métodos clásicos se ejecutan de forma independiente, llamando directamente a sus respectivos archivos:

-   `laneemdenode45.m`: resuelve la ecuación de Lane-Emden usando `ode45` de Matlab.

-   `laneemdenRK4.m`: resuelve la ecuación de Lane-Emden mediante el método de Runge-Kutta clásico de orden 4 con paso fijo.

Estos archivos no dependen de `laneemdenprueba.m` y pueden ejecutarse directamente desde la consola de Matlab proporcionando los parámetros adecuados (como el índice politrópico y el intervalo de integración).


## Ejecutar un método específico

Si se quiere ejecutar un método en particular, hay que invocarlo directamente. Por ejemplo, con el comando de Matlab:

```
> theta = laneemdenode45(n, epsilon, puntofinal);
```

en el que se tiene que consultar el comentario inicial dentro de la función para ver la sintaxis exacta de sus parámetros.


## Visualización de resultados

La visualización se realiza mediante el script:

```
> graficas
```

Al ejecutar `graficas.m`, el programa genera automáticamente **tres figuras**:

-   Una figura con los resultados obtenidos por los métodos **SLM** y **QLM**.

-   Una figura con las soluciones obtenidas mediante `ode45 y =RK4`.

-   Una figura con las soluciones analíticas de la ecuación de Lane-Emden.

Este archivo permite comparar visualmente el comportamiento de los distintos métodos numéricos y el de las soluciones exactas. Un punto a tener muy en cuenta es que para poder visualizar la gráfica de las funciones **ode45** y **RK4** es necesario leer el comentario final de ambas funciones, ya que si no se tiene en cuenta solamente representará un punto, el primer cero de la ecuación. Esto es debido a que ambas funciones están programadas para que devuelvan el punto donde se consigue el primer cero, por ello, para una buena representación es necesario el vector entero, y para ello, se necesita leer el comentario al que antes se ha hecho referencia para ejecutarlo correctamente.


## Qué puede esperar el usuario

Al utilizar este repositorio, el usuario puede:

-   Resolver numéricamente la ecuación de Lane-Emden mediante métodos iterativos y métodos clásicos.

-   Comparar gráfica y analíticamente **SLM** y **QLM**.

-   Comparar gráficamente `ode45` y `RK4`.

-   Contrastar las soluciones analíticas con las soluciones numéricas.
