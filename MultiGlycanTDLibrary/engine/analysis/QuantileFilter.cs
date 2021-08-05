using MultiGlycanTDLibrary.engine.search;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MultiGlycanTDLibrary.engine.analysis
{
    public class QuantileFilter : IFilter
    {
        double quantile = 0.25;
        double cut_off = 0;
        List<SearchResult> targets = new List<SearchResult>();
        List<SearchResult> decoys = new List<SearchResult>();
        public QuantileFilter(double quantile)
        {
            this.quantile = quantile;
        }

        public List<SearchResult> Filter()
        {
            return targets.Where(p => p.Score > cut_off).ToList();
        }

        // https://stackoverflow.com/questions/8137391/percentile-calculation
        public static double Percentile(IEnumerable<double> data, double percentile)
        {
            var elements = data.ToArray();
            Array.Sort(elements);
            double realIndex = percentile * (elements.Length - 1);
            int index = (int)realIndex;
            double frac = realIndex - index;
            if (index + 1 < elements.Length)
                return elements[index] * (1 - frac) + elements[index + 1] * frac;
            else
                return elements[index];
        }

        public void Init()
        {
            if (decoys.Count > 0)
                cut_off = Percentile(
                    decoys.Where(p => p.Score > 0).Select(p => p.Score), quantile);

        }

        public void set_data(List<SearchResult> targets, List<SearchResult> decoys)
        {
            this.targets = targets;
            this.decoys = decoys;
        }
    }
}
