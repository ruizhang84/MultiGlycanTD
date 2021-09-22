using SpectrumProcess.algorithm;
using System.Windows;

namespace MultiGlycanTD
{
    /// <summary>
    /// Interaction logic for ConfigureWindow.xaml
    /// </summary>
    public partial class ConfigureWindow : Window
    {
        public ConfigureWindow()
        {
            InitializeComponent();
            InitWindow();
        }
        public void InitWindow()
        {
            MS1Tol.Text = SearchingParameters.Access.MS1Tolerance.ToString();
            MSMS2Tol.Text = SearchingParameters.Access.MSMSTolerance.ToString();
            if (SearchingParameters.Access.MS1ToleranceBy
                == ToleranceBy.Dalton)
            {
                MS1TolByDalton.IsChecked = true;
            }
            if (SearchingParameters.Access.MS2ToleranceBy
                == ToleranceBy.Dalton)
            {
                MS2TolByDalton.IsChecked = true;
            }

            foreach (double ion in SearchingParameters.Access.Ions)
            {
                if (ion == MultiGlycanTDLibrary.util.mass.Spectrum.Proton)
                {
                    Proton.IsChecked = true;
                }
                else if (ion == MultiGlycanTDLibrary.util.mass.Spectrum.Sodium)
                {
                    Sodium.IsChecked = true;
                }
            }

            ThreadNums.Text = SearchingParameters.Access.ThreadNums.ToString();
            MaxCharge.Text = SearchingParameters.Access.MaxCharge.ToString();

            FDR.Text = (SearchingParameters.Access.FDR * 100.0).ToString();
            Coverage.Text = (SearchingParameters.Access.Coverage * 100).ToString();
            Similarity.Text = SearchingParameters.Access.Similarity.ToString();
            BinWidth.Text = SearchingParameters.Access.BinWidth.ToString();
        }

        private void Save_Click(object sender, RoutedEventArgs e)
        {
            if (SaveChanges())
            {
                SearchingParameters.Access.Update();
                Close();
            }
        }

        private void Cancel_Click(object sender, RoutedEventArgs e)
        {
            Close();
        }

        private bool SaveChanges()
        {
            return SaveTolerance() &&
                    SaveOutput() &&
                    SaveSearch();
        }

        private bool SaveSearch()
        {
            if (int.TryParse(ThreadNums.Text, out int nums))
            {
                ConfigureParameters.Access.ThreadNums = nums;
            }
            else
            {
                MessageBox.Show("Thread value is invalid!");
                return false;
            }
            if (int.TryParse(MaxCharge.Text, out int charge))
            {
                ConfigureParameters.Access.MaxCharge = charge;
            }
            else
            {
                MessageBox.Show("Max charge value is invalid!");
                return false;
            }
            if (double.TryParse(Similarity.Text, out double sim) && sim >= 0 && sim <= 1.0)
            {
                ConfigureParameters.Access.Similarity = sim;
            }
            else
            {
                MessageBox.Show("Cosine Similarity is invalid!");
                return false;
            }
            if (double.TryParse(BinWidth.Text, out double binWidth) && binWidth > 0)
            {
                ConfigureParameters.Access.BinWidth = binWidth;
            }
            else
            {
                MessageBox.Show("Spectrum BinWidth is invalid!");
                return false;
            }
            if (Proton.IsChecked == true || Sodium.IsChecked == true)
            {
                ConfigureParameters.Access.Ions.Clear();
                if (Proton.IsChecked == true)
                    ConfigureParameters.Access.Ions
                        .Add(MultiGlycanTDLibrary.util.mass.Spectrum.Proton);
                if (Sodium.IsChecked == true)
                    ConfigureParameters.Access.Ions
                        .Add(MultiGlycanTDLibrary.util.mass.Spectrum.Sodium);
            }
            else
            {
                MessageBox.Show("At least one Ion type!");
                return false;
            }

            return true;
        }


        private bool SaveTolerance()
        {
            if (double.TryParse(MS1Tol.Text, out double tol))
            {
                ConfigureParameters.Access.MS1Tolerance = tol;
            }
            else
            {
                MessageBox.Show("MS tolerance value is invalid!");
                return false;
            }

            if (double.TryParse(MSMS2Tol.Text, out tol))
            {
                ConfigureParameters.Access.MSMSTolerance = tol;
            }
            else
            {
                MessageBox.Show("MSMS tolerance value is invalid!");
                return false;
            }
            return true;
        }

        private bool SaveOutput()
        {

            if (double.TryParse(FDR.Text, out double rate) && rate >= 0 && rate <= 100)
            {
                ConfigureParameters.Access.FDR = rate * 0.01;
            }
            else
            {
                MessageBox.Show("Filter level is invalid!");
                return false;
            }
            if (double.TryParse(Coverage.Text, out double portion) && portion >= 0 && portion <= 100)
            {
                ConfigureParameters.Access.Coverage = portion * 0.01;
            }
            else
            {
                MessageBox.Show("The coverage is invalid!");
                return false;
            }
            return true;
        }

        private void MS1TolByPPM_Checked(object sender, RoutedEventArgs e)
        {
            ConfigureParameters.Access.MS1ToleranceBy = ToleranceBy.PPM;
        }

        private void MS1TolByDalton_Checked(object sender, RoutedEventArgs e)
        {
            ConfigureParameters.Access.MS1ToleranceBy = ToleranceBy.Dalton;
        }

        private void MS2TolByPPM_Checked(object sender, RoutedEventArgs e)
        {
            ConfigureParameters.Access.MS2ToleranceBy = ToleranceBy.PPM;
        }

        private void MS2TolByDalton_Checked(object sender, RoutedEventArgs e)
        {
            ConfigureParameters.Access.MS2ToleranceBy = ToleranceBy.Dalton;
        }
    }
}
