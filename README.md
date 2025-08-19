# Trabajo Práctico Final - TCyM - FIUBA
El objetivo de este trabajo práctico es desarrollar un modelo numérico en Python para simular el comportamiento térmico de un muro opaco de un edificio cuyo interior se encuentra climatizado. Se busca evaluar la carga térmica instantánea que atraviesa el cerramiento a lo largo de un ciclo diario de 24 horas, considerando la variación de las condiciones externas (temperatura del aire exterior e insolación incidente). El desarrollo teórico junto con los resultados y demás detalles se pueden encontrar en el siguiente [informe](https://github.com/agustinj-fiuba/tcym_tpf_modelado_muro/blob/main/informe.pdf).

Para este trabajo se utilizó la librería [fipy](https://www.ctcms.nist.gov/fipy/) junto con [customtkinter](https://customtkinter.tomschimansky.com/) para hacer una interfaz gráfica de usuario (GUI) de fácil acceso para la carga de datos. Además se permite la carga de archivos de tipo .json para hacer comparaciones entre diferentes muros de manera rápida. 

## Instalación

Se debe instalar [Python](https://www.python.org/downloads/). Una vez instalado, se deben instalar las librerias que figuran en *requirements.txt*, para ello se utiliza el siguiente comando:

```
pip install -r requirements.txt
```

Por último, se debe clonar el repositorio sin modificar el nombre de los archivos.
