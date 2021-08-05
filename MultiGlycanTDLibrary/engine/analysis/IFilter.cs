using MultiGlycanTDLibrary.engine.search;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MultiGlycanTDLibrary.engine.analysis
{
    public interface IFilter
    {
        void Init();
        void set_data(List<SearchResult> targets, List<SearchResult> decoys);
        List<SearchResult> Filter();
    }
}
