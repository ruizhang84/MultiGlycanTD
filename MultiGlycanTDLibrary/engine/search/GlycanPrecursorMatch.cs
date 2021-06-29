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
        protected double cutoff_ = 0.05;
        protected int top_ = 3;

        public GlycanPrecursorMatch(ISearch<string> searcher, CompdJson compdJson)
        {
            searcher_ = searcher;
            mass_map_ = compdJson.MassMap;
            distr_map_ = compdJson.DistrMap;

            List<Point<string>> glycans_ = new List<Point<string>>();
            foreach (string compose in mass_map_.Keys)
            {
                // take top 3 highest peaks with at least 0.05 distr
                var sorted = distr_map_[compose]
                    .Select((x, i) => new KeyValuePair<double, int>(x, i))
                    .OrderByDescending(x => x.Key)
                    .Take(top_).Where(x => x.Key > cutoff_)
                    .ToList();
                List<int> idx = sorted.Select(x => x.Value).ToList();

                foreach (int i in idx)
                {
                    double mass = mass_map_[compose][i];
                    Point<string> glycan = new Point<string>(mass, compose);
                    glycans_.Add(glycan);
                }
            }
            searcher_.Init(glycans_);
        }

        public List<string> Match(double precursor, int charge, double ion = 1.0078)
        {
            double mass = util.mass.Spectrum.To.Compute(precursor,ion, charge);
            return searcher_.SearchContent(mass).Distinct().ToList();
        }
    }
}
