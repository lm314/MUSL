function [ fi ] = AbsorptivePotential( q2,element, E, temperature)
%ABSORPTIVEPOTENTIAL Only includes vaues for Si, and uses the
%parameterization from Peng, Ren, Dudarev, and Whelen (1996) "Debye-Waller
%Factors and Absorptive Scattering Factors of Elemental Crystal." The
%parameterization assumes the beam energy is 100keV, so we must include a
%multiplicative factor to account for this. The input E is in Volts and q2
%is in angstroms^-2
switch element
    case 14
        switch temperature
            case -1
                fi = 0;
                return
            case 65
                a = [-0.0028 0.0119 0.0104 0.0064 0.0025 0.0381 0.1925 0.5310 1.5658 9.6306];
            case 70
                a = [-0.0027 0.0126 0.0107 0.0058 0.0023 0.0377 0.2005 0.5806 1.7611 10.5806];
            case 75
                a = [-0.0028 0.0127 0.0108 0.0058 0.0024 0.0382 0.2025 0.5845 1.7728 10.6593];
            case 77
                a = [-0.0028 0.0129 0.0108 0.0056 0.0023 0.0381 0.2058 0.6035 1.8379 10.8728];
            case 80
                a = [-0.0028 0.0132 0.0111 0.0054 0.0022 0.0383 0.2087 0.6238 1.9607 11.7245];
            case 90
                a = [-0.0029 0.0133 0.0111 0.0057 0.0023 0.0396 0.2118 0.6138 1.8803 11.1625];
            case 100
                a = [-0.0029 0.0136 0.0113 0.0058 0.0024 0.0407 0.2181 0.6275 1.9000 11.1498];
            case 110
                a = [-0.0031 0.0134 0.0115 0.0064 0.0026 0.0428 0.2187 0.5980 1.7645 10.4942];
            case 120
                a = [-0.0031 0.0149 0.0122 0.0053 0.0022 0.0429 0.2370 0.7079 2.3007 13.0434];
            case 150
                a = [-0.0032 0.0182 0.0132 0.0036 0.0020 0.0457 0.2815 0.9440 3.5756 13.4231];
            case 200
                a = [-0.0039 0.0204 0.0146 0.0051 0.0021 0.0570 0.3259 0.9670 3.4619 15.0027];
            case 250
                a = [-0.0044 0.0269 0.0156 0.0030 0.0026 0.0678 0.4239 1.3762 5.5655 12.0619];
            case 275
                a = [-0.0049 0.0268 0.0133 0.0070 0.0040 0.0745 0.4463 1.1144 2.5980 11.5706];
            case 285
                a = [-0.0050 0.0281 0.0124 0.0073 0.0044 0.0766 0.4655 1.1564 2.4371 11.0350];
            case 290
                a = [-0.0050 0.0285 0.0131 0.0070 0.0042 0.0777 0.4729 1.1851 2.6479 11.6520];
            case 293
                a = [-0.0050 0.0288 0.0108 0.0087 0.0047 0.0783 0.4792 1.1336 2.2353 10.7654];
            case 295
                a = [-0.0050 0.0298 0.0110 0.0076 0.0048 0.0785 0.4873 1.2383 2.2678 10.6100];
            case 325
                a = [-0.0001 -0.0058 0.0337 0.0175 0.0057 0.0004 0.0971 0.5278 1.7340 10.4062];
            case 400
                a = [-0.0002 -0.0070 0.0388 0.0189 0.0067 0.0012 0.1211 0.6202 1.9082 10.4631];
            case 500
                a = [-0.0070 0.0494 0.0170 -0.0052 0.0106 0.1198 0.8251 2.9204 7.9357 11.0553];
            otherwise
               errorTemp();
        end
    otherwise
        errorElem();
end
    fi = a(1).*exp(-a(6)*q2) + a(2).*exp(-a(7)*q2) + a(3).*exp(-a(8)*q2) +...
            a(4).*exp(-a(9)*q2) + a(5).*exp(-a(10)*q2);
    
    B100 = 1+1.9569341*100000;
    BE = 1+1.9569341*E;
    
    fi = B100/BE*fi;
end

function errorTemp()
     error('Error Temperature AbsorptivePotential: Parameterization does not exist at sepecfied temperature for absorptive potential.')
end

function errorElem()
    error('Error Element AbsorptivePotential: element specified is not  parameterized.')
end

