################################################################################

# Creo que esta parte del análisis (lo primero que hice, no lo que hice justo 
# antes de irme del laboratorio, que fue lo que tenía que ver con la creación
# de un clasificador que usaba los valores de las proteínas de suero) se tiene
# que dividir en 3 partes:

# a) Procesamiento de los datos de las proteínas del proteoma, y no del
# fosfoproteoma, para darme cuenta de qué datos preprocesados usé finalmente,
# ver qué figuras hice, y terminar esta sección con la creación de un
# clasificador que usaba unas 30 proteínas que se detectaban en suero sanguíneo
# por encima de un umbral de concentración mínima, que fueron las proteínas
# putativas para detectarse en las muestras de suero sanguíneo de los pacientes,
# con lo que prededir la severidad futura del covid (esto sería la segunda
# parte).

# b) Clustering de las proteínas según cómo cambiaban sus niveles de
# concentración a lo largo del tiempo con Mfuzz.

# IMPORTANTE: Puede que, para sacar las proteínas del clasificador que elaboré
# antes de irme del laboratorio, lo que hiciera fuera obtenerlas a partir de los
# datos preprocesados con R, y que no usara los datos procesados de Proteome
# Discoverer, o que sí usara los datos de Proteome Discoverer en ambos casos.
# Creo que, lo que sucedió en realidad, fue que usé los datos de R para una
# cosa, y los de Proteome Discoverer, para la otra. Pero en el artículo no
# aparece indicado así, sino que solo usé los de Proteome Discoverer

# c) Hacer lo mismo que en la parte b), pero con el fosfoproteoma.

################################################################################

## 1) Proteins clustering considering their abundance change over time



## 2) Generating a machine learning classifier based on the proteins 

# ?
# a) Loading previously generated RDA containing original protein data.
# b) Loading Proteome Discoverer normalized data.
# c) Loading Proteome Discoverer t-tests results.
# d) Loading paxdb human serum proteins concentration.
# e) Creating tags for every class from the data.
# ?

