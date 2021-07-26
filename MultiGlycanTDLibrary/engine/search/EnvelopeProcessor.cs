using SpectrumData;
using SpectrumProcess.algorithm;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MultiGlycanTDLibrary.engine.search
{
    public class EnvelopeProcessor
    {
        double range = 1; // 1 mz
        ISearch<IPeak> searcher;

        public EnvelopeProcessor(ToleranceBy by = ToleranceBy.Dalton, double tol = 0.01,
            double searchRange = 1.0)
        {
            searcher = new BucketSearch<IPeak>(by, tol);
            range = searchRange;
        }

        public void Init(List<IPeak> peaks)
        {
            searcher.Init(peaks
                .Select(p => new Point<IPeak>(p.GetMZ(), p)).ToList());
        }

        public void SetTolerance(double tol)
        {
            searcher.SetTolerance(tol);
        }

        public void SetToleranceBy(ToleranceBy by)
        {
            searcher.SetToleranceBy(by);
        }

        public SortedDictionary<int, List<IPeak>> Cluster(double mz, int charge)
        {
            //int: diff of isotope
            SortedDictionary<int, List<IPeak>> cluster =
                new SortedDictionary<int, List<IPeak>>();

            double steps = 1.0 / charge;
            int index = 0;
            while (steps * index < range) // search with 1 mz
            {
                double target = mz + steps * index;
                List<IPeak> isotopics = searcher.SearchContent(target);
                if (isotopics.Count > 0)
                    cluster[index] = isotopics;
                index++;
            }
            index = -1;
            while (steps * index > -range)
            {
                double target = mz + steps * index;
                List<IPeak> isotopics = searcher.SearchContent(target);
                if (isotopics.Count > 0)
                    cluster[index] = isotopics;
                index--;
            }

            return cluster;
        }
    }
}
