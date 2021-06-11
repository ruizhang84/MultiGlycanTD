using MultiGlycanTDLibrary.engine.analysis;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Shapes;
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
                Analyze(file, search.Target(), search.Decoy());
            }

            UpdateSignal("Done");
            return Task.CompletedTask;
        }

        private void Analyze(string msPath, List<SearchResult> targets, List<SearchResult> decoys)
        {
            //QuantileFilter filter = new QuantileFilter(SearchingParameters.Access.Quantile);
            string targetpath = System.IO.Path.Combine(System.IO.Path.GetDirectoryName(msPath),
                System.IO.Path.GetFileNameWithoutExtension(msPath) + ".csv");
            MultiThreadingSearchHelper.ReportResults(targetpath, targets);
            string decoyPath = System.IO.Path.Combine(System.IO.Path.GetDirectoryName(msPath),
                System.IO.Path.GetFileNameWithoutExtension(msPath) + "_decoy.csv");
            MultiThreadingSearchHelper.ReportResults(decoyPath, decoys);
            //FMMFDRFilter filter = new FMMFDRFilter(0.01);
            FDRFilter filter = new FDRFilter(0.01);
            filter.set_data(targets, decoys);
            filter.Init();
            List<SearchResult> results = filter.Filter();
            string path = System.IO.Path.Combine(System.IO.Path.GetDirectoryName(msPath),
                System.IO.Path.GetFileNameWithoutExtension(msPath) + "_filtered.csv");
            MultiThreadingSearchHelper.ReportResults(path, results);
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
