# Tests de regression

On n'est pas encore au point où les tests de régression sont automatiques. Ça viendra, mais je ne sais pas comment faire à l'heure actuelle.  
On pourra toujours tester les choses suivantes :

## Solvatation du méthane, mmax=0, L=30^3, N=128^3

`input.dft.in`  reads
```
boxnod = 128 128 128
boxlen = 30. 30. 30.
mmax = 0
bulk_density = 0.0333
temperature = 300.
```

```
methane, asthagiri et al.
1
#    charge     sigma     epsilon    lambda1_mol   lambda2_mol   x            y           z      Zatomic 
1    0.0     3.73      1.23        0.0           0.0       0.0         0.0          0.0      16
```

ΔG MDFT = 40.93 kJ/mol  
P-scheme correction =       -0.00 kJ/mol  
PBC correction      =       -0.00 kJ/mol  

The oxygen - methane rdf, i.e., `output/rdf.out`:
```
  0.00000000       0.00000000    
  5.85937500E-02   0.00000000    
  0.175781250       0.00000000    
  0.292968750       0.00000000    
  0.410156250       0.00000000    
  0.527343750       0.00000000    
  0.644531250       0.00000000    
  0.761718750       0.00000000    
  0.878906250       0.00000000    
  0.996093750       0.00000000    
   1.11328125       0.00000000    
   1.23046875       0.00000000    
   1.34765625       0.00000000    
   1.46484375       0.00000000    
   1.58203125       0.00000000    
   1.69921875       0.00000000    
   1.81640625       0.00000000    
   1.93359375       0.00000000    
   2.05078125       0.00000000    
   2.16796875       0.00000000    
   2.28515625       0.00000000    
   2.40234375       0.00000000    
   2.51953125       0.00000000    
   2.63671875       0.00000000    
   2.75390625       0.00000000    
   2.87109375       0.00000000    
   2.98828125      0.200137332    
   3.10546875      0.990553617    
   3.22265625       1.86837518    
   3.33984375       2.16407776    
   3.45703125       2.05543756    
   3.57421875       1.80082786    
   3.69140625       1.50076342    
   3.80859375       1.25746310    
   3.92578125       1.09979069    
   4.04296875      0.986889422    
   4.16015625      0.915402889    
   4.27734375      0.872444391    
   4.39453125      0.849646747    
   4.51171875      0.841537952    
   4.62890625      0.845047355    
   4.74609375      0.856178164    
   4.86328125      0.873611331    
   4.98046875      0.896426558    
   5.09765625      0.922896743    
   5.21484375      0.950939596    
   5.33203125      0.981483579    
   5.44921875       1.01356661    
   5.56640625       1.04826605    
   5.68359375       1.08105838    
   5.80078125       1.10328507    
   5.91796875       1.10634363    
   6.03515625       1.08892012    
   6.15234375       1.05662513    
   6.26953125       1.02239776    
   6.38671875      0.993812740    
   6.50390625      0.972700477    
   6.62109375      0.960209370    
   6.73828125      0.954916179    
   6.85546875      0.955317318    
   6.97265625      0.959755301    
   7.08984375      0.966549575    
   7.20703125      0.974425256    
   7.32421875      0.983208299    
   7.44140625      0.991785526    
   7.55859375      0.998880982    
   7.67578125       1.00477183    
   7.79296875       1.00907028    
   7.91015625       1.01199555    
   8.02734375       1.01342273    
   8.14453125       1.01348865    
   8.26171875       1.01252985    
   8.37890625       1.01045156    
   8.49609375       1.00725448    
   8.61328125       1.00340962    
   8.73046875      0.999325514    
   8.84765625      0.995505393    
   8.96484375      0.992905438    
   9.08203125      0.991280854    
   9.19921875      0.991138339    
   9.31640625      0.991678357    
   9.43359375      0.993411541    
   9.55078125      0.995067298    
   9.66796875      0.997352242    
   9.78515625      0.998946190    
   9.90234375       1.00073624    
   10.0195312       1.00164866    
   10.1367188       1.00257671    
   10.2539062       1.00264370    
   10.3710938       1.00269485    
   10.4882812       1.00231481    
   10.6054688       1.00185525    
   10.7226562       1.00126135    
   10.8398438       1.00054657    
   10.9570312       1.00012219    
   11.0742188      0.999530435    
   11.1914062      0.999357283    
   11.3085938      0.998943567    
   11.4257812      0.999091506    
   11.5429688      0.999054193    
   11.6601562      0.999388516    
   11.7773438      0.999516249    
   11.8945312      0.999967098    
   12.0117188       1.00016224    
   12.1289062       1.00049376    
   12.2460938       1.00051296    
   12.3632812       1.00064576    
   12.4804688       1.00047600    
   12.5976562       1.00034273    
   12.7148438       1.00013423    
   12.8320312      0.999839306    
   12.9492188      0.999692261    
   13.0664062      0.999339938    
   13.1835938      0.999231696    
   13.3007812      0.999015749    
   13.4179688      0.999000669    
   13.5351562      0.998918772    
   13.6523438      0.999054849    
   13.7695312      0.999062896    
   13.8867188      0.999288976    
   14.0039062      0.999438584    
   14.1210938      0.999659359    
   14.2382812      0.999763548    
   14.3554688      0.999998927    
   14.4726562       1.00003111    
   14.5898438       1.00016654    
   14.7070312       1.00016236    
   14.8242188       1.00019217    
   14.9414062       1.00009918    
```