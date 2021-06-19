using NUnit.Framework;
using SpectrumData;
using SpectrumData.Reader;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace NUnitTestProject
{
    public class SpectrumSimUnitTest
    {
        public double computeDot(List<IPeak> l1, List<IPeak> l2)
        {
            if (l1.Count == 0 || l2.Count == 0)
            {
                return 0;
            }
            int count = Math.Min(l1.Count, l2.Count);
            List<IPeak> t1 = l1.OrderByDescending(x => x.GetIntensity())
                .Take(count).OrderBy(x => x.GetMZ()).ToList();
            List<IPeak> t2 = l2.OrderByDescending(x => x.GetIntensity()).
                Take(count).OrderBy(x => x.GetMZ()).ToList();
            double numerator = 0;
            for (int i = 0; i < t1.Count; i++)
            {
                numerator += t1[i].GetIntensity() * t2[i].GetIntensity();
            }
            return numerator;
        }

        public double computeCos(List<IPeak> p1, List<IPeak> p2, double tol)
        {
            double lowerBound = Math.Min(p1.Min(x => x.GetMZ()), p2.Min(x => x.GetMZ()));
            double upperBound = Math.Max(p1.Max(x => x.GetMZ()), p2.Max(x => x.GetMZ()));
            int bucketNums = (int)Math.Ceiling((upperBound - lowerBound + 1) / tol);

            List<IPeak>[] q1 = new List<IPeak>[bucketNums];
            List<IPeak>[] q2 = new List<IPeak>[bucketNums];

            for (int i = 0; i < bucketNums; i++)
            {
                q1[i] = new List<IPeak>();
                q2[i] = new List<IPeak>();
            }

            foreach (IPeak pk in p1)
            {
                int index = (int)Math.Ceiling((pk.GetMZ() - lowerBound) / tol);
                q1[index].Add(pk);
            }
            foreach (IPeak pk in p2)
            {
                int index = (int)Math.Ceiling((pk.GetMZ() - lowerBound) / tol);
                q2[index].Add(pk);
            }

            double numerator = 0;
            for (int i = 0; i < bucketNums; i++)
            {
                numerator += computeDot(q1[i], q2[i]);
            }

            double denominator1 = 0;
            foreach (IPeak pk in p1)
            {
                denominator1 += pk.GetIntensity() * pk.GetIntensity();
            }
            double denominator2 = 0;
            foreach (IPeak pk in p2)
            {
                denominator2 += pk.GetIntensity() * pk.GetIntensity();
            }
            double denominator = Math.Sqrt(denominator1) * Math.Sqrt(denominator2);
            return numerator / denominator;
        }
        [Test]
        public void CosTest()
        {
            // read spectrum
            string path = @"C:\Users\iruiz\Downloads\MSMS\HBS1_dextrinspkd_C18_10252018.raw";
            ThermoRawSpectrumReader reader = new ThermoRawSpectrumReader();
            reader.Init(path);

            ISpectrum A = reader.GetSpectrum(2657);
            ISpectrum B = reader.GetSpectrum(3438);

            Console.WriteLine(computeCos(A.GetPeaks(), B.GetPeaks(), 0.1));

            Assert.Pass();
        }
    }
}
