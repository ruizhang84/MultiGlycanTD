using MultiGlycanTDLibrary.algorithm;
using MultiGlycanTDLibrary.model;
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
        // distributions
        Dictionary<string, List<double>> mass_map_;
        Dictionary<string, List<double>> distr_map_;
        protected double cutoff_;

        public GlycanPrecursorMatch(ISearch<string> searcher, CompdJson compdJson,
            double cutoff = 0.01)
        {
            searcher_ = searcher;
            distr_map_ = compdJson.DistrMap;
            mass_map_ = compdJson.MassMap;
            cutoff_ = cutoff;

            List<Point<string>> glycans_ = new List<Point<string>>();
            foreach (string compose in mass_map_.Keys)
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

        public List<string> Match(double precursor, int charge, double ion = 1.0078)
        {
            double mass = util.mass.Spectrum.To.Compute(precursor,ion, charge);
            return searcher_.Search(mass).Distinct().ToList();
        }
    }
}
