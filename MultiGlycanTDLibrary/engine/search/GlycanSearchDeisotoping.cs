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
    using GlycanFragments = Dictionary<FragmentTypes, List<string>>;

    public class GlycanSearchDeisotoping : IGlycanSearch
    {
        ISearch<GlycanFragments> searcher_;
        AveragineDeisotoping deisotoping_;
        Dictionary<string, List<string>> id_map_;
        readonly int maxCharge = 3; // it is not likely a higher charge for fragments.
        readonly int minMatches = 5; // it is not likely only match a few peaks.

        public GlycanSearchDeisotoping(
            ISearch<GlycanFragments> searcher, 
            GlycanJson glycanJson,
            AveragineDeisotoping deisotoping)
        {
            searcher_ = searcher;
            List<Point<GlycanFragments>> points
                = new List<Point<GlycanFragments>>();
            foreach (double mass in glycanJson.FragmentMap.Keys)
            {
                points.Add(new Point<GlycanFragments>
                    (mass, glycanJson.FragmentMap[mass]));
            }
            searcher_.Init(points);
            id_map_ = glycanJson.IDMap;
            deisotoping_ = deisotoping;
        }

        protected void UpdateMatch(PeakMatch match, IPeak peaks,
            FragmentTypes type, int potentialMatches, double expectMZ)
        {
            double diff = Math.Abs(expectMZ - peaks.GetMZ());
            double prevDiff = Math.Abs(match.TheoreticMZ - peaks.GetMZ());
            if (prevDiff > diff)
            {
                match.TheoreticMZ = expectMZ;
                match.Potentials = potentialMatches;
                match.IonTypes.Clear();
                match.IonTypes.Add(type);
            }
            else if (prevDiff == diff)
            {
                match.Potentials = Math.Min(match.Potentials, potentialMatches);
                match.IonTypes.Add(type);
            }

        }

        protected List<SearchResult> PickTop(
            Dictionary<string, SearchResult> results)
        {
            List<SearchResult> topResults = new List<SearchResult>();
            foreach (string glycan in results.Keys)
            {
                SearchResult result = results[glycan];
                if (result.Matches.Count < minMatches)
                    continue;
                topResults.Add(result);
            }

            return topResults;
        }

        protected void SearchPeaks(int index, List<DeisotopingPeak> peaks,
            double ion, int charge,
            Dictionary<string, string> glycanCandid, 
            Dictionary<string, SearchResult> results)
        {
            IPeak peak = peaks[index];
            double mass = util.mass.Spectrum.To.Compute(peak.GetMZ(),
                       ion, charge);
            List<Point<GlycanFragments>> glycans = searcher_.Search(mass);

            // make records
            foreach (Point<GlycanFragments> pt in glycans)
            {
                GlycanFragments fragments = pt.Content();
                double expectMZ = util.mass.Spectrum.To.ComputeMZ(pt.Value(), ion, charge);
                foreach (FragmentTypes type in fragments.Keys)
                {
                    foreach (string glycan in fragments[type])
                    {
                        if (!glycanCandid.ContainsKey(glycan))
                        {
                            continue;
                        }

                        if (!results.ContainsKey(glycan))
                        {
                            results[glycan] = new SearchResult();
                            results[glycan].Ion = ion;
                            results[glycan].Glycan = glycan;
                            results[glycan].Composition = glycanCandid[glycan];
                        }

                        if (!results[glycan].Matches.ContainsKey(index))
                        {
                            results[glycan].Matches[index] = new PeakMatch();
                            results[glycan].Matches[index].Peak = peaks[index];
                        }

                        UpdateMatch(results[glycan].Matches[index],
                            peaks[index], type, fragments[type].Count, expectMZ);
                    }
                }
            }
        }

        public List<SearchResult> Search(List<string> candidates, List<IPeak> peaks,
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
            List<DeisotopingPeak> deisotopingPeaks =  deisotoping_.Process(peaks, ion);

            // search peaks glycan_id->peak_index
            Dictionary<string, SearchResult> results
                = new Dictionary<string, SearchResult>();
            for (int i = 0; i < deisotopingPeaks.Count; i++)
            {
                DeisotopingPeak peak = deisotopingPeaks[i];
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