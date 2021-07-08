using MultiGlycanTDLibrary.engine.score;
using NUnit.Framework;
using SpectrumData;
using SpectrumData.Reader;
using SpectrumProcess.algorithm;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace NUnitTestProject
{
    public class SpectrumClusterUnitTest
    {
        [Test]
        public void CusterPeaksTest()
        {
            // read spectrum
            string path = @"C:\Users\iruiz\Downloads\MSMS\122123_13_C18_120min_60oC_50cm.raw";
            string output = @"C:\Users\iruiz\Downloads\MSMS\peak_clusters.csv";
            ThermoRawSpectrumReader reader = new ThermoRawSpectrumReader();
            reader.Init(path);

            int scan = 17450;
            int k = 4;
            int maxIter = 1000;
            double tol = 0.01;
            ClusterKMeans<IPeak> cluster = new ClusterKMeans<IPeak>(k, maxIter, tol);
            Dictionary<double, int> counts = new Dictionary<double, int>();
            

            using (FileStream ostrm = new FileStream(output, FileMode.OpenOrCreate, FileAccess.Write))
            {
                using (StreamWriter writer = new StreamWriter(ostrm))
                {
                    writer.WriteLine("mz,intensity,cluster,min");
                    int i = scan;

                    List<IPeak> peaks = reader.GetSpectrum(i).GetPeaks();
                    List<Point<IPeak>> points =
                    peaks.Select(p => new Point<IPeak>(p.GetIntensity(), p)).ToList();
                    cluster.Run(points);
                    double minClusterIntensity = int.MaxValue;
                    int minClusterIndex = 0;
                    foreach (int index in cluster.Clusters.Keys)
                    {
                        double average =
                            cluster.Clusters[index].Average(peaks => peaks.Content().GetIntensity());
                        if (average < minClusterIntensity)
                        {
                            minClusterIntensity = average;
                            minClusterIndex = index;
                        }
                    }

                    for (int index = 0; index < peaks.Count; index++)
                    {
                        string outputString = peaks[index].GetMZ().ToString() + "," +
                            peaks[index].GetIntensity().ToString() + ","; 
                                    
                        int clusterIndex = cluster.Index[index];
                        outputString += clusterIndex.ToString() + ",";
                        if (clusterIndex == minClusterIndex)
                        {
                            outputString += "1";
                        }
                        else
                        {
                            outputString += "0";
                        }
                        writer.WriteLine(outputString);
                    }
                        

                    
                }
            }
        }
    }
}
