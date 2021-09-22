using MultiGlycanTDLibrary.engine.score;
using MultiGlycanTDLibrary.engine.search;
using SpectrumData;
using SpectrumProcess.algorithm;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MultiGlycanTDLibrary.engine.annotation
{
    public class CoverageValidator
    {

        GlycanAnnotationSearcher searcher_;
        protected ClusterKMeans<IPeak> cluster;
        double cut_off = 0.5;

        public CoverageValidator(GlycanAnnotationSearcher searcher, 
            double cut_off, int k = 3, int maxIter = 1000, double tol = 0.01)
        {
            searcher_ = searcher;
            this.cut_off = cut_off;
            cluster = new ClusterKMeans<IPeak>(k, maxIter, tol);
        }

        public bool Valid(List<IPeak> peaks, SearchResult result)
        {
            // init
            searcher_.InitAnnotation(result);
            // cluster major peaks
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

            // search peaks
            int matched = 0;
            int nTotal = 0;
            for (int index = 0; index < peaks.Count; index++)
            {
                int clusterIndex = cluster.Index[index];
                if (clusterIndex == minClusterIndex)
                    continue;

                IPeak peak = peaks[index];
                for (int charge = 1; charge <= result.Charge; charge++)
                {
                    double mass = util.mass.Spectrum.To.Compute(peak.GetMZ(), result.Ion, charge);
                    List<GlycanAnnotated> glycans = searcher_.Search(mass);
                    if (glycans.Count > 0)
                    {
                        matched++;
                        break;
                    }
                }
                nTotal++;
            }

            double coverage = matched * 1.0 / nTotal;
            return coverage > cut_off;
        }
    }
}
