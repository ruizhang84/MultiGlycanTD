using SpectrumData;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MultiGlycanTDLibrary.engine.search
{
    public interface IGlycanSearch
    {
        List<SearchResult> Search(List<string> candidates, List<IPeak> peaks,
            int precursorCharge, double ion);
    }
}
