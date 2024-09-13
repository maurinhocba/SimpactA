# **Guía para el Uso de GitHub, Mantenimiento de Versiones y Flujos de Trabajo**

## **1\. Uso de GitHub y Branches para Mantenimiento de Versiones**

### **1.1. Estructura de Branches**

* **`main`**: Branch principal donde se encuentran las versiones estables que están listas para producción.  
* **`develop`**: Branch secundario donde se integran las nuevas características y correcciones antes de ser incluidas en `main`.  
* **Branches de Características** (`feature/nombre-caracteristica`):  
  * Usados para desarrollar nuevas funcionalidades.  
  * Ejemplo: `feature/nueva-caracteristica`.  
  * Se crean a partir de `develop` y, al finalizar, se fusionan de nuevo en `develop`.  
* **Branches de Corrección de Errores** (`bugfix/nombre-bug`):  
  * Usados para corregir errores específicos.  
  * Ejemplo: `bugfix/correccion-error`.  
  * Se crean a partir de `develop` y, al finalizar, se fusionan de nuevo en `develop`.  
* **Branches de Lanzamiento** (`release/vX.Y.Z`):  
  * Usados para preparar una nueva versión antes de fusionarla en `main`.  
  * Ejemplo: `release/v1.0.0`.  
  * Permiten pruebas finales y corrección de errores antes de ser integradas en `main`.  
* **Branches de Hotfix** (`hotfix/nombre-hotfix`):  
  * Usados para correcciones urgentes que deben aplicarse directamente en `main`.  
  * Ejemplo: `hotfix/correccion-urgente`.  
  * Se crean a partir de `main` y luego se fusionan tanto en `main` como en `develop`.

### **1.2. Ejemplo de Comandos para Manejo de Branches**

* Crear un branch de desarrollo:

bash  
`git checkout -b develop`

* Crear un branch para una nueva característica a partir de `develop`:

bash  
`git checkout -b feature/nueva-caracteristica develop`

* Crear un branch para una corrección de errores a partir de `develop`:

bash  
`git checkout -b bugfix/correccion-error develop`

* Crear un branch para preparar una nueva versión:

bash  
`git checkout -b release/v1.0.0 develop`

* Crear un branch para una corrección urgente:

bash  
`git checkout -b hotfix/correccion-urgente main`

## **2\. Recomendaciones para Nombrar Versiones**

### **2.1. Versionado Semántico (SemVer)**

* **Formato:** `MAJOR.MINOR.PATCH`  
  * **MAJOR**: Cambios incompatibles en la API.  
  * **MINOR**: Nuevas características compatibles con versiones anteriores.  
  * **PATCH**: Correcciones de errores compatibles.  
* **Ejemplo**: `1.3.2`  
  * `1`: Cambio mayor.  
  * `3`: Tres nuevas versiones menores añadidas.  
  * `2`: Dos correcciones de errores.

### **2.2. Uso de Etiquetas en Git**

* Crear una etiqueta para una nueva versión:

bash  
`git tag -a v1.0.0 -m "Descripción de la versión 1.0.0"`

* Listar etiquetas existentes:

bash  
`git tag`

* Pujar la etiqueta al repositorio remoto:

bash  
`git push origin v1.0.0`

### **2.3. Nombres de Branches y Convenciones**

* **`develop`**: Branch para desarrollo.  
* **`feature/nueva-caracteristica`**: Branch para nuevas funcionalidades.  
* **`bugfix/correccion-error`**: Branch para correcciones de errores.  
* **`release/v1.0.0`**: Branch para preparar una nueva versión.  
* **`hotfix/correccion-urgente`**: Branch para correcciones urgentes.

### **2.4. Versiones Preliminares y de Prueba**

* **Alpha:** `1.0.0-alpha` (versiones iniciales, no estables)  
* **Beta:** `1.0.0-beta` (versiones estables, posibles errores)  
* **RC (Release Candidate):** `1.0.0-rc1` (candidatas a la versión final)

### **2.5. Documentación de Cambios**

* Mantén un archivo `CHANGELOG.md` que documente los cambios en cada versión:

markdown  
`# Changelog`

`## [1.0.0] - 2024-08-01`  
`### Added`  
`- Nueva característica X.`  
`- Mejora en el rendimiento de Y.`

`### Fixed`  
`- Corregido error en Z.`

`## [0.9.0] - 2024-07-15`  
`### Added`  
`- Primera versión beta con funcionalidades básicas.`

### **2.6. Codificación de Fechas**

Si necesitas versiones basadas en fechas:

* **Formato:** `YYYYMMDD`  
* **Ejemplo:** `1.0.0-20240824`

---

