# Entrega 1. Resumen artículo del proyecto  

**Integrantes:**  
- Contreras Cornejo Yeimi Gissel  
- Mendoza Suárez Marina  
- Velazquez Aldrete Zyanya Valentina  

---

## Información sobre el organismo de referencia  

*Klebsiella* es un género de bacterias Gram-negativas, perteneciente a la familia *Enterobacteriaceae*. Son bacilos encapsulados, no móviles, facultativos anaerobios. Son ubicuas en el ambiente (suelo, agua) y forman parte de la microbiota intestinal normal de humanos y animales.  

Actualmente, las bacterias del género *Klebsiella* causan frecuentemente infecciones nosocomiales en humanos. Los principales reservorios patógenos para la transmisión de *Klebsiella* son el tracto gastrointestinal y las manos del personal hospitalario. Debido a su rápida propagación en el entorno hospitalario, estas bacterias tienden a causar brotes nosocomiales.  *(Podschun & Ullmann, 1998)*  

Sin embargo, cepas aisladas de medios ambientales presentan un aliado potencial. Como se menciona en el artículo, *Klebsiella sp.* se utiliza como agente de remediación en la investigación para la eliminación de colorantes azoicos, fenoles, hexahidro-1,3,5-trinitro-1,3,5-triazina, insecticidas organoclorados, aminas aromáticas y muchas otras sustancias tóxicas.  

No obstante, aún no se han estudiado la biodegradación de diferentes fuentes de cianuro, la optimización del medio que contiene cianuro ni el efecto de diferentes iones en la biodegradación del cianuro utilizando *K. pneumoniae*. Ya se ha determinado la eficiencia de biodegradación de *K. pneumoniae* en presencia de diferentes compuestos de cianuro:  

- Trihidrato de hexacianoferrato(II) de potasio  
- Hidrato de tetracianonickelato(II) de potasio  
- Decahidrato de ferrocianuro de sodio  

*(Avcioglu & Bilkay, 2016)*  

La biorremediación que *Klebsiella* ejerce podría aplicarse a varios metales pesados, pero en este caso, hay uno que resulta de nuestro interés debido a su uso masivo en la producción de acero inoxidable, pigmentos o cromado industrial.  

---

## Cromo hexavalente (Cr(VI))  

El Cr(VI) es extremadamente tóxico y cancerígeno, soluble y móvil, con gran capacidad de contaminación en mantos acuíferos.  

- El ion cromato es isoeléctrico e isoestructural con los iones sulfato y fosfato (esenciales para la célula).  
- Por ello, puede acceder masivamente a las células.  
- Su alto potencial de reducción positivo inicia una reacción en cadena generando intermediarios inestables (CrV y CrIV).  
- Estos reaccionan con oxígeno molecular o peróxido de hidrógeno, produciendo una avalancha de ROS.  

La estructura final de Cr(VI) es el **cromo trivalente (Cr(III))**, el cual:  
- Tiene alta afinidad por grupos donadores de electrones (bases del DNA, carboxilos de proteínas).  
- Forma enlaces covalentes muy fuertes con el DNA, creando aductos y enlaces cruzados DNA-proteína.  
- Bloquea la replicación, impide la transcripción e interfiere en la reparación del DNA.  

Si la reducción ocurre **fuera de la célula**, ya sea por agentes químicos naturales o enzimas bacterianas (como en *Klebsiella*), el producto final es insoluble:  

- El Cr(III) precipita como **hidróxido de cromo (Cr(OH)₃)**.  
- Forma un sólido inerte que no penetra la célula.  

*Klebsiella* ha evolucionado para realizar esta reducción extracelular o en compartimentos especializados, como estrategia de **detoxificación y supervivencia**. Esto le permite:  
- Vivir en ambientes letales.  
- Evitar que los metales pesados ingresen al citoplasma y causen daño.  

---

## Objetivos  

### Principal  
- Evaluar el **perfil transcriptómico global** de células adaptadas al Cr(VI) frente a controles sin Cr(VI), centrándose en genes diferencialmente expresados implicados en respuestas al estrés y procesos metabólicos.  
- Comprender los mecanismos de **resistencia y adaptación al Cr(VI)**.  

### Secundarios  
- Analizar y verificar el **potencial de biorremediación** de la cepa *Klebsiella* frente al Cr(VI).  
- Analizar la **interacción de respuestas transcripcionales** que dotan a *Klebsiella* de capacidad reductora de Cr(VI).  
- Entender cómo esta bacteria logra **reducir Cr(VI) a Cr(III)** y tolerar concentraciones muy altas.  

---

## Resultados esperados  

A partir del análisis diferencial con datos similares al artículo *“Transcriptome Analysis Reveals Cr(VI) Adaptation Mechanisms in Klebsiella sp. Strain AqSCr”*:  

- Se espera identificar un **conjunto significativo de genes diferencialmente expresados** en respuesta a Cr(VI).  
  - Genes sobreexpresados (~255).  
  - Genes sub-expresados (~240).  

- Estos genes se agruparían en **categorías funcionales clave**, reflejando las estrategias de adaptación:  
  - **Respuesta al estrés oxidativo:**  
    - ↑ *sodA* (superóxido dismutasa de manganeso).  
    - ↑ *cybB* (superóxido oxidasa).  
    - ↓ *sodB* (superóxido dismutasa de hierro).  
  - **Captación de nutrientes esenciales:** azufre, hierro, fosfato, molibdato, manganeso.  
  - **Respuesta al estrés de la envoltura celular y osmótico:** implicación de factores sigma alternativos (*fecI, rpoE, rpoS*).  
  - **Alteraciones metabólicas:**  
    - Metabolismo de ácidos grasos.  
    - Funciones ribosomales.  
    - Reprogramación del metabolismo energético.  

Esto resalta la **complejidad de la adaptación bacteriana** frente a contaminantes como el Cr(VI).  

---

## Datos de referencia  

- **BioProject ID:** 341863  
- **Genoma de referencia:** *Klebsiella pneumoniae* NTUH-K2044  
  - 92% de identidad.  
  - 85% de cobertura.  

---

## Referencias  

- Podschun, R., & Ullmann, U. (1998). *Klebsiella spp. as Nosocomial Pathogens: Epidemiology, Taxonomy, Typing Methods, and Pathogenicity Factors.* Clinical Microbiology Reviews, 11(4), 589–603. https://doi.org/10.1128/cmr.11.4.589  
- Avcioglu, N. H., & Bilkay, I. S. (2016). *Biological Treatment of Cyanide by Using Klebsiella pneumoniae Species.* Food Technology and Biotechnology, 54(4), 450–454. https://doi.org/10.17113/ftb.54.04.16.4518  

