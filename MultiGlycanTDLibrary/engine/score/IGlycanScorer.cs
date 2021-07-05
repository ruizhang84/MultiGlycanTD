using MultiGlycanTDLibrary.engine.search;
using SpectrumData;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MultiGlycanTDLibrary.engine.score
{
    public interface IGlycanScorer
    {
        public List<SearchResult> Result();
        public void Init(ConcurrentDictionary<int, ISpectrum> spectra,
            List<SearchResult> results);
        public void Run();

    }
}
