using MultiGlycanClassLibrary.algorithm;
using MultiGlycanClassLibrary.model.glycan;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MultiGlycanClassLibrary.engine.search
{
    public class GlycanSearchCreator
    {
        double tolerance;
        ToleranceBy by;
        Dictionary<string, IGlycan> glycans_map;

        public GlycanSearchCreator(ToleranceBy by, double tolerance,
            Dictionary<string, IGlycan> glycans_map)
        {
            this.by = by;
            this.tolerance = tolerance;
            this.glycans_map = glycans_map;
        }

        public GlycanSearch Create(ProducingFragments producing)
        {
            ISearch<IGlycan> searcher = new BucketSearch<IGlycan>(by, tolerance);
            return new GlycanSearch(searcher, glycans_map, producing);
        }

    }
}
