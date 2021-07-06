using MultiGlycanTDLibrary.engine.glycan;
using MultiGlycanTDLibrary.model;
using SpectrumData;
using SpectrumProcess.algorithm;
using SpectrumProcess.deisotoping;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MultiGlycanTDLibrary.engine.search
{
    using GlycanFragments = Dictionary<FragmentType, List<string>>;

    public class GlycanSearchDeisotoping : GlycanSearch, IGlycanSearch
    {
        AveragineDeisotoping deisotoping_;

        public GlycanSearchDeisotoping(
            ISearch<GlycanFragments> searcher, GlycanJson glycanJson,
            AveragineDeisotoping deisotoping): base(searcher, glycanJson)
        {
            deisotoping_ = deisotoping;
        }

        public override List<SearchResult> Search(
            List<string> candidates, List<IPeak> peaks,
            int precursorCharge, double ion = 1.0078)
        {
            // process composition, id -> compos
            Dictionary<string, string> glycanCandid = new Dictionary<string, string>();
            foreach (string composition in candidates)
            {
                foreach (string glycan in id_map_[composition])
                {
                    glycanCandid[glycan] = composition;
                }
            }
            // deisotoping
            List<IPeak> deisotopingPeaks = deisotoping_.Process(peaks, ion);

            // search peaks glycan_id->peak_index
            Dictionary<string, SearchResult> results
                = new Dictionary<string, SearchResult>();
            for (int i = 0; i < deisotopingPeaks.Count; i++)
            {
                DeisotopingPeak peak = deisotopingPeaks[i] as DeisotopingPeak;
                if (peak.ChargeAssigned())
                {
                    SearchPeaks(i, deisotopingPeaks, ion, peak.Charge, glycanCandid, results);
                }
                else
                {
                    for (int charge = 1; charge <= Math.Min(maxCharge, precursorCharge); charge++)
                    {
                        SearchPeaks(i, deisotopingPeaks, ion, charge, glycanCandid, results);
                    }
                }
            }

            // pick the top candidates
            return PickTop(results);
        }       

    }
}