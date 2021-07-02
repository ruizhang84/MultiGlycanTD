using NUnit.Framework;
using SpectrumData;
using SpectrumData.Reader;
using SpectrumProcess;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace NUnitTestProject
{
    public class SpectrumProcessTest
    {
        [Test]
        public void TestProcess()
        {
            // read spectrum
            string path = @"C:\Users\Rui Zhang\Downloads\HBS1_dextrinspkd_C18_10252018.raw";
            string output = @"C:\Users\Rui Zhang\Downloads\peaks.csv";
            ThermoRawSpectrumReader reader = new ThermoRawSpectrumReader();
            reader.Init(path);

            int scan = 291;
            ISpectrum ms1 = reader.GetSpectrum(scan);
            IProcess picking = new LocalNeighborPicking();

            ms1 = picking.Process(ms1);
            List<IPeak> ms1Peaks = ms1.GetPeaks();

            using (FileStream ostrm = new FileStream(output, FileMode.OpenOrCreate, FileAccess.Write))
            {
                using (StreamWriter writer = new StreamWriter(ostrm))
                {
                    writer.WriteLine("mz,intensity");
                    foreach (IPeak pk in ms1Peaks)
                    {
                        writer.WriteLine(pk.GetMZ().ToString() + "," +
                            pk.GetIntensity().ToString());
                    }
                }
            }
        }
    }
}