
## Análisis de Componentes Principales (ACP) en C++

Este proyecto implementa el **Análisis de Componentes Principales (ACP) en C++** para reducir la dimensionalidad de un conjunto de datos y visualizar la información en un espacio de menor dimensión. El objetivo fue comprender las diferencias en la programación de este problema en distintos lenguajes, como **Python**, que es *fuertemente tipado y dinámico*, y **C++**, que es *débilmente tipado y estático*. 

### Muestra del Proyecto

Esta muestra del proyecto presenta el **Plano Principal**, que representa la distribución de los individuos según las dos primeras componentes principales, y el **Círculo de Correlación**, que visualiza cómo las variables originales se relacionan con estas componentes a través de sus coordenadas.
![Image](https://github.com/user-attachments/assets/4f36f96a-29ad-4a9f-9994-8dee5f086362)

### Resultados Generados
- **Matriz de estandarizada (X)**
- **Matriz de correlaciones (R)**
- **Valores propios**
- **Vectores propios (V)**
- **Matriz de componentes principales (C)**
- **Matriz de calidades de individuos (Q)**
- **Matriz de coordenadas de variables (T)**
- **Vector de inercias de los ejes (I)**
- **Matriz de calidades de variables (S)**
- **Gráficos:** Plano Principal y Círculo de Correlación.

### Librerias/Programas a Instalar
- `Python`
- `Eigen`
- `Matplotlibcpp`
- `Numpy`

### Consideraciones 
Para ejecutar correctamente el código en **Visual Studio 2022** o **Visual Studio Community**, es necesario realizar algunas modificaciones. Como se mostró en la imagen anterior, el proyecto debe estar en *Release x64*, y se deben configurar las direcciones de las librerías mencionadas.

> Al abrir Proyecto > Propiedades, dirígete a las siguientes ubicaciones y agrega los siguientes paths:
- VC++ Directories
    - Include Directories
      - Agregar path del include de Python
        - Ej: C:\Python312\include  
      - Agregar path del numpy include de Python
        - Ej: C:\Python312\Lib\site-packages\numpy\_core\include
    - Library Directories
      - Agregar path de libs de Python
        - Ej: C:\Python312\libs
- C/C++
    - Additional Include Directories
      - Agregar path del Eigen
        - Ej: C:\Users\tatig\Documents\Unitec 2025\Lenguajes\Proyecto\eigen-3.4.0  
