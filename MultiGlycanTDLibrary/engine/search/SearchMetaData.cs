using MultiGlycanTDLibrary.model.glycan;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MultiGlycanTDLibrary.engine.search
{
    public class SearchMetaData
    {
        public List<SearchResult> Commit(List<SearchResult> results,
            double mz, int charge, int scan, double retention)
        {
            foreach (SearchResult r in results)
            {
                r.Scan = scan;
                r.PrecursorMZ = mz;
                r.Retention = retention;
                r.Charge = charge;
            }
            return results;
        }
    }
}
