namespace SpectrumData
{
    public interface ISpectrumReader
    {
        void Init(string fileName);
        int GetFirstScan();
        int GetLastScan();
        int GetMSnOrder(int scanNum);
        double GetRetentionTime(int scanNum);
        ISpectrum GetSpectrum(int scanNum);
        TypeOfMSActivation GetActivation(int scanNum);
        double GetPrecursorMass(int scanNum, int msOrder);
    }
}
