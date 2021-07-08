using SpectrumData;
using System.Collections.Generic;

namespace MultiGlycanTDLibrary.engine.search
{
    public interface IGlycanSearch
    {
        List<SearchResult> Search(List<string> candidates, List<IPeak> peaks,
            int precursorCharge, double ion);
    }
}
