using Microsoft.Win32;
using MultiGlycanTDLibrary.engine.glycan;
using MultiGlycanTDLibrary.model;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.Json;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;

namespace GlycanCalculator
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        protected string fileName;
        protected int hexNAc;
        protected int hex;
        protected int fuc;
        protected int neuAc;
        protected int neuGc;
        protected bool complexInclude;
        protected bool hybridInclude;
        protected bool highMannoseInclude;
        protected int order;
        protected bool permethylated;
        protected bool reduced;
        protected int precision;
        protected List<FragmentTypes> types = new List<FragmentTypes>();

        public MainWindow()
        {
            InitializeComponent();
        }

        private void JSONFileNames_Click(object sender, RoutedEventArgs e)
        {
            SaveFileDialog savefile = new SaveFileDialog();
            savefile.FileName = "database.json";
            savefile.Filter = "Json files (*.json) | *.json";

            if (savefile.ShowDialog() == true)
            {
                fileName = savefile.FileName;
            }
        }

        private bool SearchInit()
        {
            int bound = 12;
            if (int.TryParse(HexNAc.Text, out bound) && bound >= 0)
            {
                hexNAc = bound;
            }
            else
            {
                MessageBox.Show("HexNAc value is invalid!");
                return false;
            }
            if (int.TryParse(Hex.Text, out bound) && bound >= 0)
            {
                hex = bound;
            }
            else
            {
                MessageBox.Show("HexNAc value is invalid!");
                return false;
            }
            if (int.TryParse(Fuc.Text, out bound) && bound >= 0)
            {
                fuc = bound;
            }
            else
            {
                MessageBox.Show("Fuc value is invalid!");
                return false;
            }
            if (int.TryParse(NeuAc.Text, out bound) && bound >= 0)
            {
                neuAc = bound;
            }
            else
            {
                MessageBox.Show("NeuAc value is invalid!");
                return false;
            }
            if (int.TryParse(NeuGc.Text, out bound) && bound >= 0)
            {
                neuGc = bound;
            }
            else
            {
                MessageBox.Show("NeuGc value is invalid!");
                return false;
            }
            if (int.TryParse(Order.Text, out bound) && bound >= 0)
            {
                order = bound;
            }
            else
            {
                MessageBox.Show("Num Distribution value is invalid!");
                return false;
            }
            if (int.TryParse(Precision.Text, out bound) && bound >= 0)
            {
                precision = bound;
            }
            else
            {
                MessageBox.Show("Decimal value is invalid!");
                return false;
            }
            if (ComplexNGlycan.IsChecked == false &&
                HybridNGlycan.IsChecked == false && HighMannose.IsChecked == false)
            {
                MessageBox.Show("Choose at least one NGlycan type!");
                return false;
            }
            else
            {
                complexInclude = ComplexNGlycan.IsChecked == true;
                hybridInclude = HybridNGlycan.IsChecked == true;
                highMannoseInclude = HighMannose.IsChecked == true;
            }

            permethylated = Permethylated.IsChecked == true;
            reduced = Reduced.IsChecked == true;

            types.Clear();
            if (Bions.IsChecked == true)
                types.Add(FragmentTypes.B);
            if (Cions.IsChecked == true)
                types.Add(FragmentTypes.C);
            if (Yions.IsChecked == true)
                types.Add(FragmentTypes.Y);
            if (Zions.IsChecked == true)
                types.Add(FragmentTypes.Z);
            if (BYions.IsChecked == true)
                types.Add(FragmentTypes.BY);
            if (BZions.IsChecked == true)
                types.Add(FragmentTypes.BZ);
            if (CYions.IsChecked == true)
                types.Add(FragmentTypes.CY);
            if (YYions.IsChecked == true)
                types.Add(FragmentTypes.YY);
            if (YZions.IsChecked == true)
                types.Add(FragmentTypes.YZ);
            if (ZZions.IsChecked == true)
                types.Add(FragmentTypes.ZZ); 
            if (types.Count == 0)
            {
                MessageBox.Show("Select at least one ions!");
                return false;
            }
            
            if (fileName is null || fileName.Length == 0)
            {
                MessageBox.Show("Name the saving file!");
                return false;
            }
            return true;
        }

        private Task Process()
        {
            GlycanBuilder glycanBuilder =
                new GlycanBuilder(hexNAc, hex, fuc, neuAc, neuGc,
                complexInclude, hybridInclude, highMannoseInclude, 
                order, permethylated, reduced);
            glycanBuilder.Build();

            // distribution maps
            var distr_map = glycanBuilder.GlycanDistribMaps();
            var mass_map = glycanBuilder.GlycanMassMaps();
            CompdJson compdJson = new CompdJson()
            {
                DistrMap = distr_map,
                MassMap = mass_map
            };

            // fragmentation maps
            object obj = new object();
            Dictionary<double, List<string>> fragments = new Dictionary<double, List<string>>();
            var map = glycanBuilder.GlycanMaps();

            GlycanIonsBuilder.Build.Permethylated = permethylated;
            GlycanIonsBuilder.Build.Reduced = reduced;
            GlycanIonsBuilder.Build.Types = types;
            Parallel.ForEach(map, pair =>
            {
                var id = pair.Key;
                var glycan = pair.Value;
                if (glycan.IsValid())
                {
                    List<double> massList = GlycanIonsBuilder.Build.Fragments(glycan)
                                        .OrderBy(m => m).Select(m => Math.Round(m, 4)).ToList();
                    lock (obj)
                    {
                        foreach (double mass in massList)
                        {
                            if (!fragments.ContainsKey(mass))
                                fragments[mass] = new List<string>();
                            fragments[mass].Add(id);
                        }
                    }
                }
            });
            var composition_map = glycanBuilder.GlycanCompositionMaps();
            Dictionary<string, List<string>> id_map = new Dictionary<string, List<string>>();
            foreach (var pair in composition_map)
            {
                id_map[pair.Key] = pair.Value.Select(p => p.ID()).ToList();
            }

            // serialization
            GlycanJson glycanJson = new GlycanJson()
            {
                Compound = compdJson,
                IDMap = id_map,
                Fragments = fragments
            };
            string jsonString = JsonSerializer.Serialize(glycanJson);
            File.WriteAllText(fileName, jsonString);

            return Task.CompletedTask;
        }

        private void Search_Click(object sender, RoutedEventArgs e)
        {
            if (!SearchInit())
                return;

            Search.IsEnabled = false;
            Task.Run(Process);
        }
    }
}
