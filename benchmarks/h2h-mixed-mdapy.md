# d-SEAMS CHILL+ against mdapy and OVITO

Same three frames, cutoff 3.5 A, one thread. Thresholds were not
moved. The comparison operator was not moved either: d-SEAMS already
uses the inclusive edges in the papers, and on these frames the
inclusive and strict readings give the same label for every oxygen.

## Published inequalities

Nguyen and Molinero, J. Phys. Chem. B 2015, 119, 9369, equation (1),
averages the spherical harmonic over the four closest neighbours:

    q_lm(i) = (1/4) sum_{j=1}^{4} Y_lm(r_ij)

Equation (2) is the normalised correlation of that vector,

    c(i,j) = sum_m q_lm(i) q_lm(j)*
             / sqrt( (sum_m |q_lm(i)|^2) (sum_m |q_lm(j)|^2) )

The text says q_lm(i) projects the orientational structure of the
four closest neighbours, and that c(i,j) is the dot product of q_l
between a pair of neighbouring molecules. Table 2 then sets the
CHILL+ windows with closed ends:

    eclipsed    0.25 >= c(i,j) >= -0.35
    staggered   c(i,j) <= -0.8

The same page restates CHILL, citing Moore, de la Llave, Welke,
Scherlis and Molinero, Phys. Chem. Chem. Phys. 2010, 12, 4124:
"CHILL identifies a staggered bond when c(i, j) <= -0.8 and an
eclipsed bond when -0.05 >= c(i, j) >= -0.2." Table 2 prints the
same CHILL row, with the footnote that it is from that 2010 paper.
Both inequalities are inclusive. d-SEAMS classifyBond with the
published flag uses c <= -0.8 and -0.35 <= c <= 0.25 for CHILL+,
and c <= -0.8 and -0.2 <= c <= -0.05 for CHILL. That is the paper.
The operator was left as it is.

mdapy's CHILL+ page states the strict cuts c < -0.8 and
-0.35 < c < 0.25, and it sums q_3m over the neighbours inside the
cutoff rather than over the four closest. OVITO 3.17.2 ChillPlus
is the modifier in the ovito Python module, cutoff 3.5. d-SEAMS
also has an ovito-chill+ path (open edges, q3 on the whole cutoff
row). On these frames that path, mdapy, and the OVITO modifier
agree on every oxygen. The open and closed readings of the same
c_ij also agree on every oxygen, so the split is not the operator.

## Agreement

| frame | n | d-SEAMS vs OVITO | mdapy vs OVITO | d-SEAMS ovito-rule vs OVITO |
|---|---:|---:|---:|---:|
| interface | 512 | 20 | 0 | 0 |
| stacking | 512 | 0 | 0 | 0 |
| melted | 512 | 17 | 0 | 0 |

d-SEAMS vs OVITO is 20, 0 and 17 oxygens. mdapy vs OVITO is 0 on
every frame. The d-SEAMS ovito-chill+ path also matches the
modifier on every oxygen. Stacking has no disagreement: every
oxygen there has four neighbours inside 3.5 A, so the four
closest neighbours and the cutoff shell are the same set.

## Where the 37 labels move

Each disagreeing oxygen itself has four neighbours inside 3.5 A,
and that set is the same in d-SEAMS, mdapy and OVITO. The
correlation still moves because equation (1) and the cutoff sum
disagree on a neighbour that has more than four oxygens inside
3.5 A. d-SEAMS keeps the four closest in
q_lm. OVITO and mdapy fold the extra oxygen into q_lm. c_ij then
shifts by about 1e-2 to several tenths, and a bond that was
eclipsed becomes out of range, or the other way around. The
molecule row in Table 1 follows the new staggered and eclipsed
counts. The interfacial neighbour test is not the split: turning
it off does not change a label once the q3 shell is fixed.

No c_ij of a disagreeing oxygen lands on -0.8, -0.35 or 0.25.
For each bond below, the two letters are the inclusive class
and the strict class. They are the same letter in every row.
d-SEAMS applies the inclusive form. mdapy's documentation
applies the strict form. OVITO's modifier matches either form
of the cutoff-shell correlation on these frames. `paper` is the
four-closest c_ij. `shell` is the cutoff-shell c_ij.

### interface

oxygen 71: d-SEAMS hexagonal_ice, OVITO and mdapy interfacial_ice. partners whose shell is not 4: 69:6.
    j=82 paper=-1.000000000000 S/S shell=-1.000000000000 S/S |dc|=0.000e+00 dist_paper=2.000e-01 dist_shell=2.000e-01
    j=70 paper=-1.000000000000 S/S shell=-1.000000000000 S/S |dc|=0.000e+00 dist_paper=2.000e-01 dist_shell=2.000e-01
    j=69 paper=-0.212733152541 E/E shell=-0.412280262744 O/O |dc|=1.995e-01 dist_paper=1.373e-01 dist_shell=6.228e-02 changes class
    j=6 paper=-1.000000000000 S/S shell=-1.000000000000 S/S |dc|=0.000e+00 dist_paper=2.000e-01 dist_shell=2.000e-01

oxygen 108: d-SEAMS other, OVITO and mdapy interfacial_hydrate. partners whose shell is not 4: 288:6, 381:7.
    j=288 paper=-0.810426984660 S/S shell=-0.703055612180 O/O |dc|=1.074e-01 dist_paper=1.043e-02 dist_shell=9.694e-02 changes class
    j=367 paper=-0.181109714571 E/E shell=-0.181109714571 E/E |dc|=0.000e+00 dist_paper=1.689e-01 dist_shell=1.689e-01
    j=106 paper=0.005126573436 E/E shell=0.005126573436 E/E |dc|=0.000e+00 dist_paper=2.449e-01 dist_shell=2.449e-01
    j=381 paper=-0.410336430221 O/O shell=-0.255239541760 E/E |dc|=1.551e-01 dist_paper=6.034e-02 dist_shell=9.476e-02 changes class

oxygen 122: d-SEAMS interfacial_ice, OVITO and mdapy hexagonal_ice. partners whose shell is not 4: 124:7.
    j=123 paper=-1.000000000000 S/S shell=-1.000000000000 S/S |dc|=0.000e+00 dist_paper=2.000e-01 dist_shell=2.000e-01
    j=124 paper=-0.580438148071 O/O shell=-0.228233720810 E/E |dc|=3.522e-01 dist_paper=2.196e-01 dist_shell=1.218e-01 changes class
    j=187 paper=-1.000000000000 S/S shell=-1.000000000000 S/S |dc|=0.000e+00 dist_paper=2.000e-01 dist_shell=2.000e-01
    j=111 paper=-1.000000000000 S/S shell=-1.000000000000 S/S |dc|=0.000e+00 dist_paper=2.000e-01 dist_shell=2.000e-01

oxygen 151: d-SEAMS hexagonal_ice, OVITO and mdapy interfacial_ice. partners whose shell is not 4: 149:6.
    j=162 paper=-1.000000000000 S/S shell=-1.000000000000 S/S |dc|=0.000e+00 dist_paper=2.000e-01 dist_shell=2.000e-01
    j=150 paper=-1.000000000000 S/S shell=-1.000000000000 S/S |dc|=0.000e+00 dist_paper=2.000e-01 dist_shell=2.000e-01
    j=149 paper=-0.324068519184 E/E shell=-0.448447217564 O/O |dc|=1.244e-01 dist_paper=2.593e-02 dist_shell=9.845e-02 changes class
    j=86 paper=-1.000000000000 S/S shell=-1.000000000000 S/S |dc|=0.000e+00 dist_paper=2.000e-01 dist_shell=2.000e-01

oxygen 174: d-SEAMS hexagonal_ice, OVITO and mdapy interfacial_ice. partners whose shell is not 4: 176:5.
    j=175 paper=-1.000000000000 S/S shell=-1.000000000000 S/S |dc|=0.000e+00 dist_paper=2.000e-01 dist_shell=2.000e-01
    j=176 paper=-0.155206967309 E/E shell=-0.395319493391 O/O |dc|=2.401e-01 dist_paper=1.948e-01 dist_shell=4.532e-02 changes class
    j=239 paper=-1.000000000000 S/S shell=-1.000000000000 S/S |dc|=0.000e+00 dist_paper=2.000e-01 dist_shell=2.000e-01
    j=235 paper=-1.000000000000 S/S shell=-1.000000000000 S/S |dc|=0.000e+00 dist_paper=2.000e-01 dist_shell=2.000e-01

oxygen 225: d-SEAMS interfacial_hydrate, OVITO and mdapy other. partners whose shell is not 4: 417:6, 407:11, 304:5.
    j=417 paper=0.072933720897 E/E shell=-0.196450911646 E/E |dc|=2.694e-01 dist_paper=1.771e-01 dist_shell=1.535e-01
    j=227 paper=-0.560449928917 O/O shell=-0.560449928917 O/O |dc|=0.000e+00 dist_paper=2.104e-01 dist_shell=2.104e-01
    j=407 paper=-0.312676268025 E/E shell=-0.403908691072 O/O |dc|=9.123e-02 dist_paper=3.732e-02 dist_shell=5.391e-02 changes class
    j=304 paper=-0.206481796273 E/E shell=-0.547439809467 O/O |dc|=3.410e-01 dist_paper=1.435e-01 dist_shell=1.974e-01 changes class

oxygen 266: d-SEAMS interfacial_hydrate, OVITO and mdapy hydrate. partners whose shell is not 4: 344:6, 296:8, 448:7, 290:6.
    j=344 paper=-0.085194697212 E/E shell=-0.349144659499 E/E |dc|=2.639e-01 dist_paper=2.648e-01 dist_shell=8.553e-04
    j=296 paper=-0.123418452064 E/E shell=-0.300585698803 E/E |dc|=1.772e-01 dist_paper=2.266e-01 dist_shell=4.941e-02
    j=448 paper=-0.611578519735 O/O shell=-0.342096832435 E/E |dc|=2.695e-01 dist_paper=1.884e-01 dist_shell=7.903e-03 changes class
    j=290 paper=-0.037415510041 E/E shell=-0.046291028691 E/E |dc|=8.876e-03 dist_paper=2.874e-01 dist_shell=2.963e-01

oxygen 281: d-SEAMS other, OVITO and mdapy interfacial_hydrate. partners whose shell is not 4: 404:5, 447:8, 319:7.
    j=404 paper=0.265111921300 O/O shell=0.022214641977 E/E |dc|=2.429e-01 dist_paper=1.511e-02 dist_shell=2.278e-01 changes class
    j=393 paper=-0.135627943625 E/E shell=-0.135627943625 E/E |dc|=0.000e+00 dist_paper=2.144e-01 dist_shell=2.144e-01
    j=447 paper=-0.679265531059 O/O shell=-0.651235094603 O/O |dc|=2.803e-02 dist_paper=1.207e-01 dist_shell=1.488e-01
    j=319 paper=0.095245885344 E/E shell=-0.239722648312 E/E |dc|=3.350e-01 dist_paper=1.548e-01 dist_shell=1.103e-01

oxygen 308: d-SEAMS other, OVITO and mdapy interfacial_hydrate. partners whose shell is not 4: 387:7, 439:5, 415:5.
    j=387 paper=-0.070767111670 E/E shell=-0.014138651494 E/E |dc|=5.663e-02 dist_paper=2.792e-01 dist_shell=2.641e-01
    j=439 paper=-0.132448120425 E/E shell=-0.340267248420 E/E |dc|=2.078e-01 dist_paper=2.176e-01 dist_shell=9.733e-03
    j=415 paper=-0.652761932305 O/O shell=-0.268417801967 E/E |dc|=3.843e-01 dist_paper=1.472e-01 dist_shell=8.158e-02 changes class
    j=278 paper=-0.467877948339 O/O shell=-0.467877948339 O/O |dc|=0.000e+00 dist_paper=1.179e-01 dist_shell=1.179e-01

oxygen 320: d-SEAMS interfacial_hydrate, OVITO and mdapy other. partners whose shell is not 4: 476:3, 270:5, 283:6, 358:10.
    j=476 paper=-0.076043741052 E/E shell=-0.076043741052 E/E |dc|=0.000e+00 dist_paper=2.740e-01 dist_shell=2.740e-01
    j=270 paper=-0.304999975714 E/E shell=-0.489454825884 O/O |dc|=1.845e-01 dist_paper=4.500e-02 dist_shell=1.395e-01 changes class
    j=283 paper=-0.143695902619 E/E shell=-0.769650850781 O/O |dc|=6.260e-01 dist_paper=2.063e-01 dist_shell=3.035e-02 changes class
    j=358 paper=-0.516324889801 O/O shell=-0.431657176444 O/O |dc|=8.467e-02 dist_paper=1.663e-01 dist_shell=8.166e-02

oxygen 353: d-SEAMS other, OVITO and mdapy interfacial_hydrate. partners whose shell is not 4: 485:7, 329:5.
    j=1 paper=-0.481508254050 O/O shell=-0.481508254050 O/O |dc|=0.000e+00 dist_paper=1.315e-01 dist_shell=1.315e-01
    j=5 paper=-0.286173099578 E/E shell=-0.286173099578 E/E |dc|=0.000e+00 dist_paper=6.383e-02 dist_shell=6.383e-02
    j=485 paper=-0.438327994248 O/O shell=-0.271623695939 E/E |dc|=1.667e-01 dist_paper=8.833e-02 dist_shell=7.838e-02 changes class
    j=329 paper=-0.078389715048 E/E shell=-0.347613778168 E/E |dc|=2.692e-01 dist_paper=2.716e-01 dist_shell=2.386e-03

oxygen 365: d-SEAMS hydrate, OVITO and mdapy other. partners whose shell is not 4: 419:8, 285:10, 472:8.
    j=386 paper=0.166409949234 E/E shell=0.166409949234 E/E |dc|=0.000e+00 dist_paper=8.359e-02 dist_shell=8.359e-02
    j=419 paper=-0.324590833924 E/E shell=-0.324770005933 E/E |dc|=1.792e-04 dist_paper=2.541e-02 dist_shell=2.523e-02
    j=285 paper=-0.325833768351 E/E shell=-0.353268832558 O/O |dc|=2.744e-02 dist_paper=2.417e-02 dist_shell=3.269e-03 changes class
    j=472 paper=-0.313503238118 E/E shell=-0.651510822640 O/O |dc|=3.380e-01 dist_paper=3.650e-02 dist_shell=1.485e-01 changes class

oxygen 386: d-SEAMS interfacial_hydrate, OVITO and mdapy other. partners whose shell is not 4: 419:8, 472:8, 285:10.
    j=365 paper=0.166409949234 E/E shell=0.166409949234 E/E |dc|=0.000e+00 dist_paper=8.359e-02 dist_shell=8.359e-02
    j=419 paper=0.257874185107 O/O shell=0.289394095086 O/O |dc|=3.152e-02 dist_paper=7.874e-03 dist_shell=3.939e-02
    j=472 paper=0.068525717686 E/E shell=-0.378695396872 O/O |dc|=4.472e-01 dist_paper=1.815e-01 dist_shell=2.870e-02 changes class
    j=285 paper=-0.328292438561 E/E shell=-0.628554705518 O/O |dc|=3.003e-01 dist_paper=2.171e-02 dist_shell=1.714e-01 changes class

oxygen 393: d-SEAMS interfacial_hydrate, OVITO and mdapy hydrate. partners whose shell is not 4: 404:5, 319:7.
    j=281 paper=-0.135627943625 E/E shell=-0.135627943625 E/E |dc|=0.000e+00 dist_paper=2.144e-01 dist_shell=2.144e-01
    j=404 paper=-0.477308154419 O/O shell=-0.329979605008 E/E |dc|=1.473e-01 dist_paper=1.273e-01 dist_shell=2.002e-02 changes class
    j=502 paper=-0.236416316038 E/E shell=-0.236416316038 E/E |dc|=0.000e+00 dist_paper=1.136e-01 dist_shell=1.136e-01
    j=319 paper=-0.109780601002 E/E shell=-0.051564211753 E/E |dc|=5.822e-02 dist_paper=2.402e-01 dist_shell=2.984e-01

oxygen 413: d-SEAMS other, OVITO and mdapy interfacial_hydrate. partners whose shell is not 4: 443:5, 402:7, 395:9.
    j=443 paper=-0.500338386938 O/O shell=-0.415094472261 O/O |dc|=8.524e-02 dist_paper=1.503e-01 dist_shell=6.509e-02
    j=492 paper=-0.151897985937 E/E shell=-0.151897985937 E/E |dc|=0.000e+00 dist_paper=1.981e-01 dist_shell=1.981e-01
    j=402 paper=0.627235108125 O/O shell=-0.154716060735 E/E |dc|=7.820e-01 dist_paper=3.772e-01 dist_shell=1.953e-01 changes class
    j=395 paper=0.170455080596 E/E shell=-0.038836273449 E/E |dc|=2.093e-01 dist_paper=7.954e-02 dist_shell=2.888e-01

oxygen 435: d-SEAMS interfacial_hydrate, OVITO and mdapy hydrate. partners whose shell is not 4: 292:3, 398:7, 315:9, 483:5.
    j=292 paper=-0.143916757220 E/E shell=-0.143916757220 E/E |dc|=0.000e+00 dist_paper=2.061e-01 dist_shell=2.061e-01
    j=398 paper=0.241569778458 E/E shell=-0.123799197383 E/E |dc|=3.654e-01 dist_paper=8.430e-03 dist_shell=2.262e-01
    j=315 paper=0.265655525350 O/O shell=0.206205243731 E/E |dc|=5.945e-02 dist_paper=1.566e-02 dist_shell=4.379e-02 changes class
    j=483 paper=-0.152937717861 E/E shell=0.069607230592 E/E |dc|=2.225e-01 dist_paper=1.971e-01 dist_shell=1.804e-01

oxygen 452: d-SEAMS other, OVITO and mdapy hydrate. partners whose shell is not 4: 260:3, 443:5, 310:8, 447:8.
    j=260 paper=-0.220925690067 E/E shell=-0.220925690067 E/E |dc|=0.000e+00 dist_paper=1.291e-01 dist_shell=1.291e-01
    j=443 paper=0.131071866576 E/E shell=0.178193779943 E/E |dc|=4.712e-02 dist_paper=1.189e-01 dist_shell=7.181e-02
    j=310 paper=0.478870024899 O/O shell=0.029337752340 E/E |dc|=4.495e-01 dist_paper=2.289e-01 dist_shell=2.207e-01 changes class
    j=447 paper=0.499708206697 O/O shell=0.028373872399 E/E |dc|=4.713e-01 dist_paper=2.497e-01 dist_shell=2.216e-01 changes class

oxygen 492: d-SEAMS other, OVITO and mdapy interfacial_hydrate. partners whose shell is not 4: 429:6, 392:7, 318:5.
    j=413 paper=-0.151897985937 E/E shell=-0.151897985937 E/E |dc|=0.000e+00 dist_paper=1.981e-01 dist_shell=1.981e-01
    j=429 paper=-0.350501987271 O/O shell=-0.184081819644 E/E |dc|=1.664e-01 dist_paper=5.020e-04 dist_shell=1.659e-01 changes class
    j=392 paper=-0.166300518222 E/E shell=-0.497045378761 O/O |dc|=3.307e-01 dist_paper=1.837e-01 dist_shell=1.470e-01 changes class
    j=318 paper=0.257765239942 O/O shell=0.024517888592 E/E |dc|=2.332e-01 dist_paper=7.765e-03 dist_shell=2.255e-01 changes class

oxygen 499: d-SEAMS interfacial_hydrate, OVITO and mdapy hydrate. partners whose shell is not 4: 408:6, 463:6, 419:8.
    j=408 paper=-0.352540999053 O/O shell=0.148094963073 E/E |dc|=5.006e-01 dist_paper=2.541e-03 dist_shell=1.019e-01 changes class
    j=463 paper=-0.249748832513 E/E shell=0.177158844113 E/E |dc|=4.269e-01 dist_paper=1.003e-01 dist_shell=7.284e-02
    j=419 paper=-0.141741359151 E/E shell=-0.136204444022 E/E |dc|=5.537e-03 dist_paper=2.083e-01 dist_shell=2.138e-01
    j=396 paper=0.054845812115 E/E shell=0.054845812115 E/E |dc|=0.000e+00 dist_paper=1.952e-01 dist_shell=1.952e-01

oxygen 512: d-SEAMS interfacial_hydrate, OVITO and mdapy other. partners whose shell is not 4: 336:6.
    j=374 paper=-0.289672902535 E/E shell=-0.289672902535 E/E |dc|=0.000e+00 dist_paper=6.033e-02 dist_shell=6.033e-02
    j=12 paper=0.135385762239 E/E shell=0.135385762239 E/E |dc|=0.000e+00 dist_paper=1.146e-01 dist_shell=1.146e-01
    j=336 paper=-0.066097243109 E/E shell=-0.530477523169 O/O |dc|=4.644e-01 dist_paper=2.839e-01 dist_shell=1.805e-01 changes class
    j=324 paper=-0.488775583684 O/O shell=-0.488775583684 O/O |dc|=0.000e+00 dist_paper=1.388e-01 dist_shell=1.388e-01

### stacking

No disagreeing oxygen.

### melted

oxygen 48: d-SEAMS interfacial_hydrate, OVITO and mdapy other. partners whose shell is not 4: 465:5, 460:7.
    j=465 paper=0.181540605474 E/E shell=0.366717189065 O/O |dc|=1.852e-01 dist_paper=6.846e-02 dist_shell=1.167e-01 changes class
    j=467 paper=-0.716380859852 O/O shell=-0.716380859852 O/O |dc|=0.000e+00 dist_paper=8.362e-02 dist_shell=8.362e-02
    j=179 paper=-0.254465282442 E/E shell=-0.254465282442 E/E |dc|=0.000e+00 dist_paper=9.553e-02 dist_shell=9.553e-02
    j=460 paper=-0.178394642664 E/E shell=-0.340919226502 E/E |dc|=1.625e-01 dist_paper=1.716e-01 dist_shell=9.081e-03

oxygen 52: d-SEAMS hydrate, OVITO and mdapy interfacial_hydrate. partners whose shell is not 4: 174:6, 181:3.
    j=174 paper=-0.338789975693 E/E shell=-0.694879884549 O/O |dc|=3.561e-01 dist_paper=1.121e-02 dist_shell=1.051e-01 changes class
    j=181 paper=0.127671133412 E/E shell=0.127671133412 E/E |dc|=0.000e+00 dist_paper=1.223e-01 dist_shell=1.223e-01
    j=297 paper=0.162648489669 E/E shell=0.162648489669 E/E |dc|=0.000e+00 dist_paper=8.735e-02 dist_shell=8.735e-02
    j=169 paper=-0.243497355883 E/E shell=-0.243497355883 E/E |dc|=0.000e+00 dist_paper=1.065e-01 dist_shell=1.065e-01

oxygen 53: d-SEAMS other, OVITO and mdapy interfacial_hydrate. partners whose shell is not 4: 61:8, 182:7, 176:7.
    j=304 paper=0.238659975922 E/E shell=0.238659975922 E/E |dc|=0.000e+00 dist_paper=1.134e-02 dist_shell=1.134e-02
    j=61 paper=-0.385824185713 O/O shell=-0.112182485237 E/E |dc|=2.736e-01 dist_paper=3.582e-02 dist_shell=2.378e-01 changes class
    j=182 paper=0.367959637255 O/O shell=0.448664781272 O/O |dc|=8.071e-02 dist_paper=1.180e-01 dist_shell=1.987e-01
    j=176 paper=0.001754313414 E/E shell=0.039169096175 E/E |dc|=3.741e-02 dist_paper=2.482e-01 dist_shell=2.108e-01

oxygen 81: d-SEAMS interfacial_ice, OVITO and mdapy other. partners whose shell is not 4: 83:5, 460:7, 178:10.
    j=83 paper=-0.657284520617 O/O shell=-0.756829090978 O/O |dc|=9.954e-02 dist_paper=1.427e-01 dist_shell=4.317e-02
    j=460 paper=-0.872980583728 S/S shell=-0.797051791393 O/O |dc|=7.593e-02 dist_paper=7.298e-02 dist_shell=2.948e-03 changes class
    j=464 paper=-0.872980583728 S/S shell=-0.872980583728 S/S |dc|=0.000e+00 dist_paper=7.298e-02 dist_shell=7.298e-02
    j=178 paper=0.854275401202 O/O shell=0.420948442187 O/O |dc|=4.333e-01 dist_paper=6.043e-01 dist_shell=1.709e-01

oxygen 171: d-SEAMS interfacial_hydrate, OVITO and mdapy other. partners whose shell is not 4: 209:5, 337:5.
    j=298 paper=0.492549047121 O/O shell=0.492549047121 O/O |dc|=0.000e+00 dist_paper=2.425e-01 dist_shell=2.425e-01
    j=212 paper=-0.235138805270 E/E shell=-0.235138805270 E/E |dc|=0.000e+00 dist_paper=1.149e-01 dist_shell=1.149e-01
    j=209 paper=-0.266995714524 E/E shell=-0.444073999651 O/O |dc|=1.771e-01 dist_paper=8.300e-02 dist_shell=9.407e-02 changes class
    j=337 paper=0.033739538698 E/E shell=-0.051645093196 E/E |dc|=8.538e-02 dist_paper=2.163e-01 dist_shell=2.984e-01

oxygen 186: d-SEAMS hexagonal_ice, OVITO and mdapy interfacial_ice. partners whose shell is not 4: 315:7.
    j=188 paper=-0.111111111111 E/E shell=-0.111111111111 E/E |dc|=0.000e+00 dist_paper=2.389e-01 dist_shell=2.389e-01
    j=159 paper=-1.000000000000 S/S shell=-1.000000000000 S/S |dc|=0.000e+00 dist_paper=2.000e-01 dist_shell=2.000e-01
    j=187 paper=-0.861494563954 S/S shell=-0.861494563954 S/S |dc|=0.000e+00 dist_paper=6.149e-02 dist_shell=6.149e-02
    j=315 paper=-0.871944427620 S/S shell=-0.551053068573 O/O |dc|=3.209e-01 dist_paper=7.194e-02 dist_shell=2.011e-01 changes class

oxygen 206: d-SEAMS interfacial_ice, OVITO and mdapy hexagonal_ice. partners whose shell is not 4: 208:5.
    j=207 paper=-1.000000000000 S/S shell=-1.000000000000 S/S |dc|=0.000e+00 dist_paper=2.000e-01 dist_shell=2.000e-01
    j=208 paper=-0.374071772354 O/O shell=-0.131859819394 E/E |dc|=2.422e-01 dist_paper=2.407e-02 dist_shell=2.181e-01 changes class
    j=331 paper=-0.859889994189 S/S shell=-0.859889994189 S/S |dc|=0.000e+00 dist_paper=5.989e-02 dist_shell=5.989e-02
    j=335 paper=-1.000000000000 S/S shell=-1.000000000000 S/S |dc|=0.000e+00 dist_paper=2.000e-01 dist_shell=2.000e-01

oxygen 212: d-SEAMS other, OVITO and mdapy interfacial_hydrate. partners whose shell is not 4: 44:8, 209:5.
    j=298 paper=-0.020953276091 E/E shell=-0.020953276091 E/E |dc|=0.000e+00 dist_paper=2.710e-01 dist_shell=2.710e-01
    j=171 paper=-0.235138805270 E/E shell=-0.235138805270 E/E |dc|=0.000e+00 dist_paper=1.149e-01 dist_shell=1.149e-01
    j=44 paper=-0.628953497487 O/O shell=-0.817342585982 S/S |dc|=1.884e-01 dist_paper=1.710e-01 dist_shell=1.734e-02 changes class
    j=209 paper=0.402197400620 O/O shell=0.227623341130 E/E |dc|=1.746e-01 dist_paper=1.522e-01 dist_shell=2.238e-02 changes class

oxygen 272: d-SEAMS interfacial_ice, OVITO and mdapy other. partners whose shell is not 4: 405:6.
    j=270 paper=-0.158470627035 E/E shell=-0.158470627035 E/E |dc|=0.000e+00 dist_paper=1.915e-01 dist_shell=1.915e-01
    j=401 paper=-0.938217070303 S/S shell=-0.938217070303 S/S |dc|=0.000e+00 dist_paper=1.382e-01 dist_shell=1.382e-01
    j=405 paper=-0.938217070303 S/S shell=-0.795484278738 O/O |dc|=1.427e-01 dist_paper=1.382e-01 dist_shell=4.516e-03 changes class
    j=307 paper=-0.357202241692 O/O shell=-0.357202241692 O/O |dc|=0.000e+00 dist_paper=7.202e-03 dist_shell=7.202e-03

oxygen 278: d-SEAMS interfacial_ice, OVITO and mdapy other. partners whose shell is not 4: 407:5, 311:5.
    j=280 paper=0.101130856565 E/E shell=0.101130856565 E/E |dc|=0.000e+00 dist_paper=1.489e-01 dist_shell=1.489e-01
    j=403 paper=-0.874065467761 S/S shell=-0.874065467761 S/S |dc|=0.000e+00 dist_paper=7.407e-02 dist_shell=7.407e-02
    j=407 paper=-0.874065467761 S/S shell=-0.724574930503 O/O |dc|=1.495e-01 dist_paper=7.407e-02 dist_shell=7.543e-02 changes class
    j=311 paper=-0.407376905279 O/O shell=-0.495602734135 O/O |dc|=8.823e-02 dist_paper=5.738e-02 dist_shell=1.456e-01

oxygen 304: d-SEAMS interfacial_hydrate, OVITO and mdapy other. partners whose shell is not 4: 182:7, 176:7, 149:5.
    j=53 paper=0.238659975922 E/E shell=0.238659975922 E/E |dc|=0.000e+00 dist_paper=1.134e-02 dist_shell=1.134e-02
    j=182 paper=0.326665872437 O/O shell=0.390362693368 O/O |dc|=6.370e-02 dist_paper=7.667e-02 dist_shell=1.404e-01
    j=176 paper=0.208084151646 E/E shell=-0.209730499218 E/E |dc|=4.178e-01 dist_paper=4.192e-02 dist_shell=1.403e-01
    j=149 paper=-0.180362214585 E/E shell=-0.401864551041 O/O |dc|=2.215e-01 dist_paper=1.696e-01 dist_shell=5.186e-02 changes class

oxygen 314: d-SEAMS hexagonal_ice, OVITO and mdapy interfacial_ice. partners whose shell is not 4: 315:7.
    j=316 paper=-0.111111111111 E/E shell=-0.111111111111 E/E |dc|=0.000e+00 dist_paper=2.389e-01 dist_shell=2.389e-01
    j=287 paper=-1.000000000000 S/S shell=-1.000000000000 S/S |dc|=0.000e+00 dist_paper=2.000e-01 dist_shell=2.000e-01
    j=315 paper=-0.871944427620 S/S shell=-0.551053068573 O/O |dc|=3.209e-01 dist_paper=7.194e-02 dist_shell=2.011e-01 changes class
    j=443 paper=-1.000000000000 S/S shell=-1.000000000000 S/S |dc|=0.000e+00 dist_paper=2.000e-01 dist_shell=2.000e-01

oxygen 336: d-SEAMS hexagonal_ice, OVITO and mdapy interfacial_ice. partners whose shell is not 4: 341:5, 465:5.
    j=341 paper=-1.000000000000 S/S shell=-0.946368813875 S/S |dc|=5.363e-02 dist_paper=2.000e-01 dist_shell=1.464e-01
    j=334 paper=-0.111111111111 E/E shell=-0.111111111111 E/E |dc|=0.000e+00 dist_paper=2.389e-01 dist_shell=2.389e-01
    j=469 paper=-1.000000000000 S/S shell=-1.000000000000 S/S |dc|=0.000e+00 dist_paper=2.000e-01 dist_shell=2.000e-01
    j=465 paper=-0.853179325799 S/S shell=-0.763366388998 O/O |dc|=8.981e-02 dist_paper=5.318e-02 dist_shell=3.663e-02 changes class

oxygen 345: d-SEAMS interfacial_hydrate, OVITO and mdapy other. partners whose shell is not 4: 49:7, 340:5, 216:5.
    j=49 paper=0.039610368436 E/E shell=0.333881857347 O/O |dc|=2.943e-01 dist_paper=2.104e-01 dist_shell=8.388e-02 changes class
    j=347 paper=0.063433983539 E/E shell=0.063433983539 E/E |dc|=0.000e+00 dist_paper=1.866e-01 dist_shell=1.866e-01
    j=340 paper=-0.329857929328 E/E shell=-0.478911625031 O/O |dc|=1.491e-01 dist_paper=2.014e-02 dist_shell=1.289e-01 changes class
    j=216 paper=-0.874180358523 S/S shell=-0.841804673308 S/S |dc|=3.238e-02 dist_paper=7.418e-02 dist_shell=4.180e-02

oxygen 431: d-SEAMS interfacial_ice, OVITO and mdapy other. partners whose shell is not 4: 430:5, 302:5.
    j=458 paper=-1.000000000000 S/S shell=-1.000000000000 S/S |dc|=0.000e+00 dist_paper=2.000e-01 dist_shell=2.000e-01
    j=430 paper=-0.250231074153 E/E shell=-0.114115876370 E/E |dc|=1.361e-01 dist_paper=9.977e-02 dist_shell=2.359e-01
    j=429 paper=-0.111111111111 E/E shell=-0.111111111111 E/E |dc|=0.000e+00 dist_paper=2.389e-01 dist_shell=2.389e-01
    j=302 paper=-0.873629566284 S/S shell=-0.659204042676 O/O |dc|=2.144e-01 dist_paper=7.363e-02 dist_shell=1.408e-01 changes class

oxygen 458: d-SEAMS hexagonal_ice, OVITO and mdapy interfacial_ice. partners whose shell is not 4: 460:7, 75:3.
    j=431 paper=-1.000000000000 S/S shell=-1.000000000000 S/S |dc|=0.000e+00 dist_paper=2.000e-01 dist_shell=2.000e-01
    j=460 paper=-0.111111111111 E/E shell=-0.356505171195 O/O |dc|=2.454e-01 dist_paper=2.389e-01 dist_shell=6.505e-03 changes class
    j=75 paper=-0.958314847500 S/S shell=-0.958314847500 S/S |dc|=0.000e+00 dist_paper=1.583e-01 dist_shell=1.583e-01
    j=459 paper=-1.000000000000 S/S shell=-1.000000000000 S/S |dc|=0.000e+00 dist_paper=2.000e-01 dist_shell=2.000e-01

oxygen 470: d-SEAMS hexagonal_ice, OVITO and mdapy interfacial_ice. partners whose shell is not 4: 83:5.
    j=472 paper=-0.111111111111 E/E shell=-0.111111111111 E/E |dc|=0.000e+00 dist_paper=2.389e-01 dist_shell=2.389e-01
    j=471 paper=-1.000000000000 S/S shell=-1.000000000000 S/S |dc|=0.000e+00 dist_paper=2.000e-01 dist_shell=2.000e-01
    j=83 paper=-0.862942118870 S/S shell=-0.597154745446 O/O |dc|=2.658e-01 dist_paper=6.294e-02 dist_shell=2.028e-01 changes class
    j=87 paper=-1.000000000000 S/S shell=-1.000000000000 S/S |dc|=0.000e+00 dist_paper=2.000e-01 dist_shell=2.000e-01

## Floating point

A double near 0.35 has unit in the last place 5.551e-17. The closest any c_ij of a disagreeing oxygen comes
to -0.8, -0.35 or 0.25 is 5.020e-04 (interface oxygen 492, paper -0.350501987271, shell -0.184081819644). The closest bond that actually changes class is the same gap, 5.020e-04 (interface oxygen 492). That gap is 9.043e+12 units in the last place, not roundoff.
Inclusive and strict comparisons classify every one of these
bonds the same way. The values are not sitting on the boundary.

