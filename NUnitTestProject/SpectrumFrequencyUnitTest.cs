using NUnit.Framework;
using SpectrumData;
using SpectrumData.Reader;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace NUnitTestProject
{
    public class SpectrumFrequencyUnitTest
    {
        [Test]
        public void CountPeaksTest()
        {
            // read spectrum
            string path = @"C:\Users\iruiz\Downloads\MSMS\HBS1_dextrinspkd_C18_10252018.raw";
            string output = @"C:\Users\iruiz\Downloads\MSMS\peaks.csv";
            ThermoRawSpectrumReader reader = new ThermoRawSpectrumReader();
            reader.Init(path);

            Dictionary<double, int> counts = new Dictionary<double, int>();
            for(int i = reader.GetFirstScan(); i <= reader.GetLastScan(); i++)
            {
                if (reader.GetMSnOrder(i) == 2)
                {
                    foreach (IPeak peak in reader.GetSpectrum(i).GetPeaks())
                    {
                        double mz = Math.Round(peak.GetMZ(), 1);
                        if (!counts.ContainsKey(mz))
                        {
                            counts[mz] = 0;
                        }
                        counts[mz] += 1;
                    }
                }

            }

            using (FileStream ostrm = new FileStream(output, FileMode.OpenOrCreate, FileAccess.Write))
            {
                using (StreamWriter writer = new StreamWriter(ostrm))
                {
                    writer.WriteLine("mz,#peaks");
                    foreach (double mz in counts.Keys)
                    {
                        writer.WriteLine(mz.ToString() + "," +
                            counts[mz].ToString());
                    }
                }
            }
        }
    }
}
