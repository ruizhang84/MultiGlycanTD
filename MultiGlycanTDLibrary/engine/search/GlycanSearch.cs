using MultiGlycanClassLibrary.algorithm;
using MultiGlycanClassLibrary.model.glycan;
using SpectrumData;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MultiGlycanClassLibrary.engine.search
{
    public delegate List<double> ProducingFragments(IGlycan glycan);

    public class GlycanSearch
    {
        ISearch<IGlycan> searcher_;
        Dictionary<string, IGlycan> glycans_map_;
        ProducingFragments producing;

        public GlycanSearch(ISearch<IGlycan> searcher,
            Dictionary<string, IGlycan> glycans_map,
            ProducingFragments producing)
        {
            searcher_ = searcher;
            glycans_map_ = glycans_map;
            this.producing = producing;
        }

        void Init(List<IGlycan> candidates)
        {
            List<Point<IGlycan>> points = new List<Point<IGlycan>>();
            foreach (IGlycan glycan in candidates)
            {
                List<double> fragments = producing(glycan);
                foreach (double mass in fragments)
                {
                    Point<IGlycan> point = new Point<IGlycan>(mass, glycan);
                    points.Add(point);
                }
            }
            searcher_.Init(points);
        }

        public Dictionary<string, HashSet<int>> Search(List<IPeak> peaks, int max_charge,
            List<IGlycan> candidates)
        {
            Init(candidates);

            Dictionary<string, HashSet<int>> searched = new Dictionary<string, HashSet<int>>();

            // search peaks and create glycanNodes
            for (int i = 0; i < peaks.Count; i++)
            {
                IPeak peak = peaks[i];
                for (int charge = 1; charge <= 2; charge++)   // modify change to max of 2
                {
                    double mass = util.mass.Spectrum.To.Compute(peak.GetMZ(),
                        util.mass.Spectrum.Proton, charge);
                    List<IGlycan> glycans = searcher_.Search(mass, mass);
                    if (glycans.Count > 0)
                    {
                        foreach (IGlycan glycan in glycans)
                        {
                            string id = glycan.ID();
                            if (!searched.ContainsKey(id))
                            {
                                searched[id] = new HashSet<int>();
                            }
                            searched[id].Add(i);
                        }
                    }
                }
            }

            Dictionary<string, HashSet<int>> results = new Dictionary<string, HashSet<int>>();

            // enumerate peak Nodes
            double best = 0;
            foreach (IGlycan glycan in candidates)
            {
                HashSet<int> matches = new HashSet<int>();
                foreach (string identified_glycan_id in searched.Keys)
                {
                    if (Satisify(identified_glycan_id, glycan))
                    {
                        matches.UnionWith(searched[identified_glycan_id]);
                    }
                }
                if (matches.Count > 0)
                    results[glycan.ID()] = matches;
            }
            return results;
        }

        bool Satisify(string identified_glycan_id, IGlycan glycan)
        {
            int[] identified_glycan_table = glycans_map_[identified_glycan_id].Table();
            int[] candidate_glycan_table = glycan.Table();
            if (identified_glycan_table.Count() != candidate_glycan_table.Count())
                return false;
            for (int i = 0; i < identified_glycan_table.Length; i++)
            {
                if (candidate_glycan_table[i] < identified_glycan_table[i])
                    return false;
            }

            // check terminal 
            if (glycans_map_[identified_glycan_id].Type() == GlycanType.NGlycanComplex)
            {
                for (int i = 0; i < 4; i++)
                {
                    if ((identified_glycan_table[12 + i] > 0 || identified_glycan_table[16 + i] > 0
                            ) && identified_glycan_table[4 + i] != candidate_glycan_table[4 + i])
                        return false;
                }

            }
            else if (glycans_map_[identified_glycan_id].Type() == GlycanType.NGlycanHybrid)
            {
                for (int i = 0; i < 2; i++)
                {
                    if ((identified_glycan_table[10 + i] > 0 || identified_glycan_table[12 + i] > 0
                            ) && identified_glycan_table[6 + i] != candidate_glycan_table[6 + i])
                        return false;
                }
            }

            return true;
        }

    }
}
