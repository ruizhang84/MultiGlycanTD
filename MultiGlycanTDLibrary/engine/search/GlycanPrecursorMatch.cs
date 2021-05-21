using MultiGlycanTDLibrary.algorithm;
using MultiGlycanTDLibrary.model.glycan;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MultiGlycanTDLibrary.engine.search
{
    public class GlycanPrecursorMatch
    {
        ISearch<IGlycan> searcher_;

        public GlycanPrecursorMatch(ISearch<IGlycan> searcher)
        {
            searcher_ = searcher;
        }

        public void Init(Dictionary<string, IGlycan> glycan_map)
        {
            List<Point<IGlycan>> glycans_ = new List<Point<IGlycan>>();
            foreach (string id in glycan_map.Keys)
            {
                if (glycan_map[id].IsValid())
                {
                    foreach(double mass in glycan_map[id].Mass())
                    {
                        Point<IGlycan> glycan = new Point<IGlycan>(mass, glycan_map[id]);
                        glycans_.Add(glycan);
                    }
                }
            }
            searcher_.Init(glycans_);
        }

        public List<IGlycan> Match(double precursor, int charge)
        {
            double mass = util.mass.Spectrum.To.Compute(precursor,
                util.mass.Spectrum.Proton, charge);
            return searcher_.Search(mass).Distinct().ToList();
        }
    }
}
