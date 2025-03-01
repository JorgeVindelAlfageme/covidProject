################################################################################

# Creo que esta parte del an�lisis (lo primero que hice, no lo que hice justo 
# antes de irme del laboratorio, que fue lo que ten�a que ver con la creaci�n
# de un clasificador que usaba los valores de las prote�nas de suero) se tiene
# que dividir en 3 partes:

# a) Procesamiento de los datos de las prote�nas del proteoma, y no del
# fosfoproteoma, para darme cuenta de qu� datos preprocesados us� finalmente,
# ver qu� figuras hice, y terminar esta secci�n con la creaci�n de un
# clasificador que usaba unas 30 prote�nas que se detectaban en suero sangu�neo
# por encima de un umbral de concentraci�n m�nima, que fueron las prote�nas
# putativas para detectarse en las muestras de suero sangu�neo de los pacientes,
# con lo que prededir la severidad futura del covid (esto ser�a la segunda
# parte).

# b) Clustering de las prote�nas seg�n c�mo cambiaban sus niveles de
# concentraci�n a lo largo del tiempo con Mfuzz.

# IMPORTANTE: Puede que, para sacar las prote�nas del clasificador que elabor�
# antes de irme del laboratorio, lo que hiciera fuera obtenerlas a partir de los
# datos preprocesados con R, y que no usara los datos procesados de Proteome
# Discoverer, o que s� usara los datos de Proteome Discoverer en ambos casos.
# Creo que, lo que sucedi� en realidad, fue que us� los datos de R para una
# cosa, y los de Proteome Discoverer, para la otra. Pero en el art�culo no
# aparece indicado as�, sino que solo us� los de Proteome Discoverer

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

