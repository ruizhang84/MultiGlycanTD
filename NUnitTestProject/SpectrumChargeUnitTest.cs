using NUnit.Framework;
using SpectrumData;
using SpectrumData.Reader;
using SpectrumProcess;
using SpectrumProcess.algorithm;
using SpectrumProcess.deisotoping;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace NUnitTestProject
{
    public class SpectrumChargeUnitTest
    {
        [Test]
        public void ChargeTest()
        {
            // read spectrum
            string path = @"C:\Users\iruiz\Downloads\MSMS\HBS1_dextrinspkd_C18_10252018.raw";
            ThermoRawSpectrumReader reader = new ThermoRawSpectrumReader();
            reader.Init(path);

            // search
            int targetScan = 1006; // 1385;
            double targetMZ = 1181.61; // 1327.62;

            double searchRange = 2.0;
            IProcess process = new WeightedAveraging(new LocalNeighborPicking());


            ISpectrum ms2 = reader.GetSpectrum(targetScan);
            ms2 = process.Process(ms2);


            ICharger charger = new Patterson();
            int charge = charger.Charge(ms2.GetPeaks(), targetMZ - searchRange, targetMZ + searchRange);
            Console.WriteLine(charge);

        }
    }
}
