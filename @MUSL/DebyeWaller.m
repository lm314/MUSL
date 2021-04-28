function DW = DebyeWaller(element,temperature)
%DEBYEWALLER returns the Debye Waller factor as calculated from the
%exoperimental density of states in Peng et al (1996) "Debye-Waller Factors
%and Absorptive Scattering Factors of Elemental Crystals". The tables with
%most of the values can be found in its supplement SUP82472.
%The value returned is in angstroms^2
switch element
    case 14
        switch temperature
            case -1
                DW = 0;
            case 65
                DW = 0.2064;
            case 70
                DW = 0.2088;
            case 75
                DW = 0.2113;
            case 77
                DW = 0.2124;
            case 80
                DW = 0.2141;
            case 90
                DW = 0.2201;
            case 100
                DW = 0.2268;
            case 110
                DW = 0.2342;
            case 120
                DW = 0.2423;
            case 150
                DW = 0.2709;
            case 200
                DW = 0.3310;
            case 250
                DW = 0.4138;
            case 275
                DW = 0.4551;
            case 285
                DW = 0.4717;
            case 290
                DW = 0.4799;
            case 293
                DW = 0.4849;
            case 295
                DW = 0.4882;
            case 325
                DW = 0.5378;
            case 400
                DW = 0.6619;
            case 500
                DW = 0.8271;
            otherwise
                errorTemp();
                
        end
   case 79
        switch temperature
            case -1
                DW = 0;
            case 65
                DW = 0.1378;
            case 70
                DW = 0.1484;
            case 75
                DW = 0.1590;
            case 77
                DW = 0.1632;
            case 80
                DW = 0.1696;
            case 90
                DW = 0.1908;
            case 100
                DW = 0.2120;
            case 110
                DW = 0.2332;
            case 120
                DW = 0.2544;
            case 150
                DW = 0.3179;
            case 200
                DW = 0.4236;
            case 250
                DW = 0.5292;
            case 290
                DW = 0.6135;
            case 293
                DW = 0.6198;
            case 295
                DW = 0.6240;
            case 320
                DW = 0.6766;
            case 400
                DW = 0.8445;
            case 500
                DW = 1.053;
            otherwise
                warningTemp();
                DW = 0;
        end
    otherwise
        warningElem();
        DW = 0;
end
end

function warningTemp()
    error('DebyeWaller Temperature: The temperature specified is not in the tabled data. Debye-Waller Factor has defaulted to 0.')
end

function warningElem()
    error('DebyeWaller Element: The element specified is not in the tabled data. Debye-Waller Factor has defaulted to 0.')
end