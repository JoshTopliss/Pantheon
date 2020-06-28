import java.lang.*;
public class NozzleDesignCalc{
    //    Vars     //

    //Input
    public static double MassFlowRate = 2602; //Kg/s
    public static double OFRatio = 2.2674;
    public static double ThroatArea = .69; //meters^2
    public static double ExitArea = 11.04; //meters^2
    public static double AmbientPressure = 1; //atm

    //Output
    public static double MolecularWeightGas;
    private static double GasConstant = 8314.4598;
    public static double ChamberPressure;
    public static double ChamberTemperature;
    public static double SpecificHeatRatio;
    private static double Gama;
    private static double ChamberTemp;
    public static double MError;
    public static double ArError;
    public static double C12H26Mass;
    public static double O2Mass;
    public static double H2OMass;
    public static double CO2Mass;
    public static double COMass;
    public static double MachExit;
    public static double ExitTemperature;
    public static double ExitVelocity;
    public static double ExitPressure;
    public static double Thrust;
    public static double Isp;

    public static void main(String[] args)
    {
        MolecularWeightGas = StochCalc(OFRatio);
        ChamberPressure = PCalc(MassFlowRate);
        SpecificHeatRatio = Gama;
        ChamberTemperature = ChamberTemp;
        MachExit = MachExitCalc(ExitArea, ThroatArea);
        ExitTemperature = ChamberTemperature/(1+(Gama-1)/2*MachExit*MachExit);
        ExitVelocity = MachExit * Math.sqrt(Gama*GasConstant/MolecularWeightGas*ExitTemperature);
        ExitPressure = ChamberPressure*(Math.pow(1+(Gama-1)/2*MachExit*MachExit, -(Gama/(Gama-1))));
        Thrust = (MassFlowRate*ExitVelocity + 101325*(ExitPressure - AmbientPressure)*ExitArea)/1000;
        Isp = ExitVelocity/9.80665;
        OutputData();
    }

    // Analyzes combustion and calculates the exauhst composition and the molecular weight of the exauhst)
    public static double StochCalc(double Ratio)
    {
        var intO2mol = OFRatio*170/32;
        var C12H26mol = Math.max(0, 1-intO2mol*2/25);
        C12H26Mass = C12H26mol * 170;
        var O2mol = Math.max(0, 12.22-(12+13/2)*1);
        O2Mass = O2mol * 32;
        var H2Omol = 13*(1-C12H26mol);
        H2OMass = H2Omol * 18;
        var CO2mol = 2*(intO2mol-O2mol) - 12*(1-C12H26mol) - H2Omol;
        CO2Mass = CO2mol * 44;
        var COmol = 12*(1-C12H26mol) - CO2mol;
        COMass = COmol * 28;
        return (C12H26Mass + O2Mass + H2OMass + CO2Mass + COMass)/(C12H26mol + O2mol + H2Omol + CO2mol + COmol);
    }

    // Bounded Secant Method Solver to calculate Chamber Pressure value when passed the Mass Flow Rate
    public static double PCalc(double M)
    {
        // x = Range Maximum, y = Range Minimum, z = Range Median
        double Pmax = 1000;
        double Pmin = 1;
        double Pest = 0; // arbitrary set value
        double Mmin = MFlowCalc(Pmin);
        double Mmax = MFlowCalc(Pmax);
        double Mest = 0;
        for (int i = 0; i<10; i++)
        {
            Pest = Pmin + ((Pmax-Pmin)/(Mmax-Mmin)*(M-Mmin));
            Mest = MFlowCalc(Pest);
            if (Mest > M)
            {
                Pmax = Pest;
                Mmax = Mest;
            }
            else
            {
                Pmin = Pest;
                Mmin = Mest;
            }
        }
        MError = Mest - M;
        return Pest;
    }

    // Calculates Mass Flow Rate (unit), Specific Heat Ratio, and Chamber Temperature (unit) when passed a pressure value (atm abs)
    public static double MFlowCalc(double P)
    {
        // Takes the natural log of the pressure
        double lnP = Math.log(P);
        // Calculates Gama Value when passed the ln of the chamber pressure **current structure is specific to dodecane**
        var GamaR1 = 0.00018885 * lnP*lnP - 0.0045309*lnP + 1.238;
        var GamaR2 = 0.0001213 * lnP*lnP  - 0.0047169*lnP + 1.2342;
        Gama = GamaR1 + (OFRatio - 2.2)/(2.4 - 2.2)*(GamaR2 - GamaR1);
        // Calculate Adiabatic Flame Temp Value when passed the ln of the chamber pressure **current structure is specific to dodecane**
        var Temp1 = -5.0031 * lnP*lnP  + 137.09*lnP + 2974.2;
        var Temp2 = -2.1907 * lnP*lnP  + 141.97*lnP + 3022.8;
        ChamberTemp = Temp1 + (OFRatio - 2.2)/(2.4 - 2.2)*(Temp2 - Temp1);
        // Calculates Mass Flow Rate
        double mflow = ThroatArea * 1000 * (P * 101.325) / Math.sqrt(GasConstant * ChamberTemp/ MolecularWeightGas)*(Math.sqrt(Gama))/Math.pow(((Gama + 1)/2),((Gama + 1)/(Gama - 1)/2));
        return mflow;
    }

    //Calculates the Mach speed at the nozzle exit
    public static double MachExitCalc(double Aexit, double Athroat)
    {
        double ARatio = Aexit/Athroat;
        double Machmax = 10;
        double Machmin = 1;
        double Machest = 0;
        double Fmax = Math.log(ARatioCalc(Machmax)/ARatio);
        double Fmin = Math.log(ARatioCalc(Machmin)/ARatio);
        double Fest = 0;
        for (int i = 0; i<10; i++)
        {
            Machest = Machmin - Fmin*((Machmax-Machmin)/(Fmax-Fmin));
            Fest = Math.log(ARatioCalc(Machest)/ARatio);
            if (Fest > 0)
            {
                Machmax = Machest;
                Fmax = Fest;
            }
            else
            {
                Machmin = Machest;
                Fmin = Fest;
            }
        }
        ArError = ARatio*Math.exp(Fest) - ARatio;
        return Machest;
    }

    //Calculates Expansion Ratio when passed an Exauhst Mach
    public static double ARatioCalc(double Mach)
    {
        return Math.pow((Gama+1)/2,-((Gama+1)/(2*(Gama-1))))*Math.pow(1+(Gama-1)/2*Mach*Mach, (Gama+1)/(2*(Gama-1)))/Mach;
    }

    //Outputs Data
    public static void OutputData()
    {
        System.out.println("Thrust: " + Thrust);
        System.out.println("Specific Impulse: " + Isp);
        System.out.println("Chamber Pressure: " + ChamberPressure);
        System.out.println("Chamber Temperature: " + ChamberTemperature);
        System.out.println("Exit Pressure: " + ExitPressure);
        System.out.println("Exit Temperature: " + ExitTemperature);
        System.out.println("Exit Velocity: " + ExitVelocity);
        System.out.println("Exit Mach: " + MachExit);
    }
}
