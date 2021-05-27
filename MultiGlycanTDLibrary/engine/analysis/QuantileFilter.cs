using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MultiGlycanTDLibrary.engine.analysis
{
    public class QuantileFilter
    {
        double quantile = 0.25;
        public QuantileFilter(double quantile)
        {
            this.quantile = quantile;
        }

        public double Quantile()
        { return quantile; }

        public void SetQuantile(double quantile)
        {
            if (quantile < 0 || quantile > 1)
                return;
            this.quantile = quantile;
        }

        public List<SearchResult> Filter(List<SearchResult> target)
        {
            int take = (int) Math.Ceiling((1 - quantile) * target.Count);

            return target.OrderByDescending(p => p.Score())
                .Take(take).OrderBy(p => p.Scan()).ToList();
        }
    }
}
