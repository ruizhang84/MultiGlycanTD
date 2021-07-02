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
    public class SpectrumFragmentDeisotoping
    {
        [Test]
        public void DeisotopingTest()
        {
            // read spectrum
            string path = @"C:\Users\iruiz\Downloads\MSMS\HBS1_dextrinspkd_C18_10252018.raw";
            string output = @"C:\Users\iruiz\Downloads\MSMS\processed_peaks.csv";
            ThermoRawSpectrumReader reader = new ThermoRawSpectrumReader();
            reader.Init(path);

            // read target spectrum
            int targetScan = 1385;

            IProcess process = new WeightedAveraging(new LocalNeighborPicking());
            ISpectrum ms2 = reader.GetSpectrum(targetScan);
            ms2 = process.Process(ms2);

            // process peaks
            Averagine averagine = new Averagine(AveragineType.PermethylatedGlycan);
            AveragineDeisotoping deisotoping = new AveragineDeisotoping(averagine,
                4, ToleranceBy.Dalton, 0.1);

            List<DeisotopingPeak> peaks = deisotoping.Process(ms2.GetPeaks(), 1.0078);
            using (FileStream ostrm = new FileStream(output, FileMode.OpenOrCreate, FileAccess.Write))
            {
                using (StreamWriter writer = new StreamWriter(ostrm))
                {
                    writer.WriteLine("mz,intensity,charge");
                    foreach (DeisotopingPeak pk in peaks)
                    {
                        writer.WriteLine(pk.GetMZ().ToString() + "," +
                            pk.GetIntensity().ToString() + "," +
                            (pk.ChargeAssigned()? pk.Charge.ToString() : "0"));
                    }
                }
            }


        }
    }
}
