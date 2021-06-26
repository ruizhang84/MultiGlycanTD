using FiniteMixtureModel.ECM;
using FiniteMixtureModel.FMM;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MultiGlycanTDLibrary.engine.analysis
{
    public class FMMFDRFilter
    {
        double fdr_;
        double cutoff_;
        List<ReportResult> target_ = new List<ReportResult>();
        List<ReportResult> decoy_ = new List<ReportResult>();
        public FMMFDRFilter(double fdr)
        {
            fdr_ = fdr;
            cutoff_ = -1;
        }

        public void Init()
        {
            // init
            cutoff_ = -1;
            if (decoy_.Count == 0 || target_.Count == 0)   //trivial case
            {
                return;
            }

            List<double> target = target_.Select(p => p.Score()).ToList();
            List<double> decoy = decoy_.Select(p => p.Score()).ToList();

            // fit gamma
            ECMAlgorithm bestEcm = null;
            double bestBIC = int.MaxValue;
            for (int component = 1; component <= 4; component++)
            {
                ECMAlgorithm ecm = new ECMAlgorithm(decoy, component);
                if (!ecm.Run())
                    break;
                if (ecm.BIC() < bestBIC)
                {
                    bestBIC = ecm.BIC();
                    bestEcm = ecm;
                }
            }

            // fit GMM
            bestBIC = int.MaxValue;
            FMMAlgorithm bestFMM = null;
            for (int component = 1; component <= 10; component++)
            {
                for (double weight = 0.1; weight <= 0.5; weight += 0.1)
                {
                    FMMAlgorithm fmm = new FMMAlgorithm(target, component, bestEcm, weight);
                    fmm.Run();
                    if (fmm.BIC() < bestBIC)
                    {
                        bestFMM = fmm;
                        bestBIC = fmm.BIC();
                    }
                }

            }

            double score = 0;
            target = target.OrderByDescending(p => p).ToList();
            for (int i = 0; i < target.Count; i++)
            {
                score += bestFMM.Probability(target[i]);
                double rate = score / (i + 1);
                if (rate < fdr_)
                {
                    cutoff_ = target[i];
                }
            }
        }

        public List<ReportResult> Filter()
        {
            return target_
                .Where(p => p.Score() >= cutoff_)
                .OrderBy(p => p.Scan()).ToList();
        }

        public void set_data(List<ReportResult> targets,
            List<ReportResult> decoys)
        {
            // acquire the best score of the scan
            Dictionary<int, double> score_map = new Dictionary<int, double>();
            targets = targets.Where(p => p.Score() > 0).ToList();
            decoys = decoys.Where(p => p.Score() > 0).ToList();
            // acquire the best score of the scan
            foreach (var it in targets)
            {
                int scan = it.Scan();
                if (!score_map.ContainsKey(scan))
                {
                    score_map[scan] = it.Score();
                }
                else
                {
                    if (score_map[scan] < it.Score())
                    {
                        score_map[scan] = it.Score();
                    }
                }
            }
            // when target and decoy are in the same spectrum,
            // the one with higher score is picked.
            foreach (var it in decoys)
            {
                int scan = it.Scan();
                if (!score_map.ContainsKey(scan))
                {
                    score_map[scan] = it.Score();
                }
                else
                {
                    if (score_map[scan] < it.Score())
                    {
                        score_map[scan] = it.Score();
                    }
                }
            }

            foreach (var it in targets)
            {
                int scan = it.Scan();
                if (score_map[scan] > it.Score())
                    continue;
                target_.Add(it);
            }
            foreach (var it in decoys)
            {
                int scan = it.Scan();
                if (score_map[scan] > it.Score())
                    continue;
                decoy_.Add(it);
            }
        }

        public double Cutoff() { return cutoff_; }

        public void set_cutoff(double cutoff) { cutoff_ = cutoff; }
    }
}
