using MultiGlycanTDLibrary.engine.analysis;
using MultiGlycanTDLibrary.engine.annotation;
using MultiGlycanTDLibrary.engine.search;
using SpectrumData;
using SpectrumData.Spectrum;
using SpectrumProcess.algorithm;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Threading;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Threading;

namespace MultiGlycanTD
{
    /// <summary>
    /// Interaction logic for SearchWindow.xaml
    /// </summary>
    public partial class SearchWindow : Window
    {
        int ReadingCounter;
        int ProgressCounter;

        public SearchWindow()
        {
            InitializeComponent();
            InitProcess();
        }

        private async void InitProcess()
        {
            await Task.Run(Process);
            Close();
        }

        private Task Process()
        {

            Counter readerCounter = new Counter();
            Counter searchCounter = new Counter();
            readerCounter.progressChange += ReadProgressChanged;
            searchCounter.progressChange += SearchProgressChanged;

            int index = 1;
            ISearch<GlycanAnnotated> searcher = new BucketSearch<GlycanAnnotated>(
                SearchingParameters.Access.MS2ToleranceBy,
                SearchingParameters.Access.MSMSTolerance);
            GlycanAnnotationLazy annotator = new GlycanAnnotationLazy(searcher,
                SearchingParameters.Access.Database.Parameters);

            foreach (string file in SearchingParameters.Access.MSMSFiles)
            {
                ReadingCounter = 0;
                ProgressCounter = 0;
                UpdateProgress(100);
                Readingprogress(100);
                UpdateSignal($"Searching... ({index++}/{SearchingParameters.Access.MSMSFiles.Count})");
                MultiThreadingSearch search =
                    new MultiThreadingSearch(file, readerCounter, searchCounter,
                       SearchingParameters.Access.Database);
                search.Run();
                UpdateSignal("Analyzing...");
                ConcurrentDictionary<int, ISpectrum> spectra = search.MSMSSpectra();
                Analyze(file, spectra, search.Target(), search.Decoy(), annotator);
            }

            UpdateSignal("Done");
            return Task.CompletedTask;
        }

        private void Analyze(string msPath,
            ConcurrentDictionary<int, ISpectrum> spectra,
            List<SearchResult> targets, List<SearchResult> decoys,
            GlycanAnnotationLazy annotator)
        {
            //QuantileFilter filter = new QuantileFilter(SearchingParameters.Access.Quantile);
            string targetpath = System.IO.Path.Combine(System.IO.Path.GetDirectoryName(msPath),
                System.IO.Path.GetFileNameWithoutExtension(msPath) + "_targets.csv");
            MultiThreadingSearchHelper.Report(targetpath, targets.Where(r => r.Score > 0).ToList());
            string decoyPath = System.IO.Path.Combine(System.IO.Path.GetDirectoryName(msPath),
                System.IO.Path.GetFileNameWithoutExtension(msPath) + "_decoys.csv");
            MultiThreadingSearchHelper.Report(decoyPath, decoys.Where(r => r.Score > 0).ToList());

            //FMMFDRFilter filter = new FMMFDRFilter(0.01);
            FDRFilter filter = new FDRFilter(SearchingParameters.Access.FDR);
            filter.set_data(targets, decoys);
            filter.Init();
            List<SearchResult> results = filter.Filter();
            string path = System.IO.Path.Combine(System.IO.Path.GetDirectoryName(msPath),
                System.IO.Path.GetFileNameWithoutExtension(msPath) + "_filtered.csv");
            MultiThreadingSearchHelper.Report(path, results);

            //Annotation
            Dictionary<int, List<PeakAnnotated>> annotations =
                new Dictionary<int, List<PeakAnnotated>>();
            foreach (SearchResult result in results)
            {
                MS2Spectrum spectrum = spectra[result.Scan] as MS2Spectrum;
                annotations[result.Scan] = new List<PeakAnnotated>();
                foreach (double ion in SearchingParameters.Access.Ions)
                {
                    annotations[result.Scan].AddRange(annotator.Annotated(spectrum.GetPeaks(), result));
                }
            }
            string annotatedPath = System.IO.Path.Combine(System.IO.Path.GetDirectoryName(msPath),
               System.IO.Path.GetFileNameWithoutExtension(msPath) + "_annotated.csv");
            MultiThreadingSearchHelper.AnnotationReport(annotatedPath, annotations);

        }


        private void UpdateSignal(string signal)
        {
            Dispatcher.BeginInvoke(
                DispatcherPriority.Normal,
                new ThreadStart(() => Signal.Text = signal));
        }

        private void UpdateProgress(int total)
        {
            Dispatcher.BeginInvoke(
                DispatcherPriority.Normal,
                new ThreadStart(() =>
                {
                    SearchingStatus.Value = ProgressCounter * 1.0 / total * 1000.0;
                }));
        }

        private void Readingprogress(int total)
        {
            Dispatcher.BeginInvoke(
                DispatcherPriority.Normal,
                new ThreadStart(() =>
                {
                    ReadingStatus.Value = ReadingCounter * 1.0 / total * 1000.0;
                }));
        }

        private void SearchProgressChanged(object sender, ProgressingEventArgs e)
        {
            Interlocked.Increment(ref ProgressCounter);
            UpdateProgress(e.Total);
        }

        private void ReadProgressChanged(object sender, ProgressingEventArgs e)
        {
            Interlocked.Increment(ref ReadingCounter);
            Readingprogress(e.Total);
        }
    }
}
