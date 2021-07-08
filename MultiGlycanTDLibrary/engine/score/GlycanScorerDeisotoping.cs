using MultiGlycanTDLibrary.engine.search;
using SpectrumData;
using SpectrumProcess.algorithm;
using SpectrumProcess.deisotoping;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace MultiGlycanTDLibrary.engine.score
{
    public class GlycanScorerDeisotoping : GlycanScorer, IGlycanScorer
    {
        Averagine Averagine;
        int MaxCharge;
        ToleranceBy By;
        double Tol;

        public GlycanScorerDeisotoping(Averagine averagine,
             int maxCharge, ToleranceBy by, double tol, int thread = 4,
             double similar = 0.9, double binWidth = 1.0)
            : base(thread, similar, binWidth)
        {
            Averagine = averagine;
            MaxCharge = maxCharge;
            By = by;
            Tol = tol;
        }

        protected void LocalAssignScore(int scan,
            AveragineDeisotoping deisotoping)
        {
            foreach (SearchResult result in SpectrumResults[scan])
            {
                List<IPeak> peaks = deisotoping
                    .Process(Spectra[scan].GetPeaks(), result.Ion)
                    .Select(p => p as IPeak).ToList();
                result.Score = GlycanScorerHelper.ComputeScore(result, peaks);
                result.Fit = GlycanScorerHelper.ComputeFit(result, peaks);
            }
        }

        protected void AssignScoreTask(ConcurrentQueue<int> ScanQueue)
        {
            AveragineDeisotoping deisotoping =
                new AveragineDeisotoping(Averagine, MaxCharge, By, Tol);
            while (ScanQueue.TryDequeue(out int scan))
            {
                LocalAssignScore(scan, deisotoping);
            }
        }

        public override void AssignScore()
        {
            ConcurrentQueue<int> ScanQueue =
                new ConcurrentQueue<int>(SpectrumResults.Keys);

            List<Task> scoer = new List<Task>();
            for (int i = 0; i < Thread; i++)
            {
                Task LastTask = new Task(() => AssignScoreTask(ScanQueue));
                LastTask.Start();
                scoer.Add(LastTask);
            }
            Task.WaitAll(scoer.ToArray());
        }

    }
}
