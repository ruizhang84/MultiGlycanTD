using NUnit.Framework;
using SpectrumData.Reader;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace NUnitTestProject
{
    public class SpectrumCountUnitTest
    {
        [Test]
        public void CountPeaksTest()
        {
            // read spectrum
            string path = @"C:\Users\Rui Zhang\Downloads\HBS1_dextrinspkd_C18_10252018.raw";
            string output = @"C:\Users\Rui Zhang\Downloads\peaks.csv";
            ThermoRawSpectrumReader reader = new ThermoRawSpectrumReader();
            reader.Init(path);

            Dictionary<int, int> counts = new Dictionary<int, int>();
            for(int i = reader.GetFirstScan(); i <= reader.GetLastScan(); i++)
            {
                if (reader.GetMSnOrder(i) == 2)
                {
                    counts[i] = reader.GetSpectrum(i).GetPeaks().Count;
                }

            }

            using (FileStream ostrm = new FileStream(output, FileMode.OpenOrCreate, FileAccess.Write))
            {
                using (StreamWriter writer = new StreamWriter(ostrm))
                {
                    writer.WriteLine("scan,#peaks");
                    foreach (int scan in counts.Keys)
                    {
                        writer.WriteLine(scan.ToString() + "," +
                            counts[scan].ToString());
                    }
                }
            }
        }
    }
}
