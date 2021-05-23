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
        ISearch<string> searcher_;
        Dictionary<string, List<IGlycan>> glycan_compound_map_;
        Dictionary<string, List<double>> mass_map_;
        Dictionary<string, List<double>> distr_map_;
        protected double cutoff_;

        public GlycanPrecursorMatch(ISearch<string> searcher,
            Dictionary<string, List<IGlycan>> glycan_compound_map,
            Dictionary<string, List<double>> mass_map,
            Dictionary<string, List<double>> distr_map,
            double cutoff = 0.01
            )
        {
            searcher_ = searcher;
            glycan_compound_map_ = glycan_compound_map;
            distr_map_ = distr_map;
            mass_map_ = mass_map;
            cutoff_ = cutoff;

            List<Point<string>> glycans_ = new List<Point<string>>();
            foreach (string compose in glycan_compound_map_.Keys)
            {
                for(int i = 0; i < mass_map_[compose].Count; i++)
                {
                    double mass = mass_map_[compose][i];
                    double distr = distr_map_[compose][i];
                    if (distr < cutoff)
                        continue;
                    Point<string> glycan = new Point<string>(mass, compose);
                    glycans_.Add(glycan);
                }
            }
            searcher_.Init(glycans_);
        }

        public List<IGlycan> Match(double precursor, int charge)
        {
            double mass = util.mass.Spectrum.To.Compute(precursor,
                util.mass.Spectrum.Proton, charge);
            return searcher_.Search(mass).Distinct()
                .SelectMany(compose => glycan_compound_map_[compose])
                .Where(g => g.IsValid())
                .ToList();
        }
    }
}
