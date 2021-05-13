using MultiGlycanClassLibrary.algorithm;
using MultiGlycanClassLibrary.model.glycan;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MultiGlycanClassLibrary.engine.search
{
    public class PrecursorMatch
    {
        ISearch<IGlycan> searcher_;
        List<Point<IGlycan>> glycans_
            = new List<Point<IGlycan>>();

        public PrecursorMatch(ISearch<IGlycan> searcher)
        {
            searcher_ = searcher;
        }

        public void Init(Dictionary<string, IGlycan> glycan_map)
        {
            foreach (var it in glycan_map.Keys)
            {
                if (glycan_map[it].IsValid())
                {
                    Point<IGlycan> glycan = new Point<IGlycan>(glycan_map[it].Mass(), glycan_map[it]);
                    glycans_.Add(glycan);
                }
            }
            searcher_.Init(glycans_);
        }

        public List<IGlycan> Match(double precursor, int charge)
        {
            double mass = util.mass.Spectrum.To.Compute(precursor,
                util.mass.Spectrum.Proton, charge);
            return searcher_.Search(mass);
        }
    }
}
