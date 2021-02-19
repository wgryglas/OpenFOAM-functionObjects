

# Turbulence Statistics tool

This is desigend for usage with DNS/LES simulation. The turbulence statistics
are gathered during runtime basing on velocity/presssure fileds, without any
turbulence model - just a raw flow data.


## Usage

file: controlDictionary
```
functionObjectName
{
  libs ("<USER COMPILATION LOCATION>/platforms/linux64GccDPInt32Opt/lib/libuserFunctionObjects.so"); // or just "libs ("libuserFunctionObjects.so");" if default user libraries location used
  type turbStatistics;
  nu 1.8e-5;

  rhoRef 1.0; //optional for compressible flow. Not required for incompressible. It will overwrite density field
  rhoName rho; //optionl for compressible flow if different variable should be used as density
}
```

file: fvSchemes
```
...

  divSchemes
  {
     ...

     statistics Gauss linear;

     ...
  }

...

```
or set "Gauss linear" for default scheme:
```
...
  divSchemes
  {
     ...

     default Gauss linear;

     ...
  }
...
```
