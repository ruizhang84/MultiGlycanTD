using System;
using System.Collections.Generic;
using System.Linq;

namespace SpectrumProcess.algorithm
{
    public class BucketSearch<T> : ISearch<T>
    {
        protected double lower_;
        protected double upper_;
        protected double tolerance_;
        protected ToleranceBy type_;
        protected List<List<Point<T>>> data_ = new List<List<Point<T>>>();

        public BucketSearch(ToleranceBy by, double tol)
        {
            tolerance_ = tol;
            type_ = by;
        }

        // init without points
        public void Init()
        {
            lower_ = 100;
            upper_ = 20000;
            if (type_ == ToleranceBy.PPM)
                PPMInit();
            else if (type_ == ToleranceBy.Dalton)
                DaltonInit();
        }

        void DaltonInit()
        {
            // fill the bucket
            int size = (int)Math.Ceiling((upper_ - lower_ + 1.0) / tolerance_);
            size = Math.Max(size, 0);

            data_ = new List<List<Point<T>>>(size);
            for (int i = 0; i < size; i++)
            {
                data_.Add(new List<Point<T>>());
            }
        }

        void PPMInit()
        {
            // fill the bucket
            double ratio = 1.0 / (1.0 - tolerance_ / 1000000);
            int size = (int)Math.Ceiling(Math.Log(upper_ / lower_) / Math.Log(ratio));
            size = Math.Max(size, 0);
            data_ = new List<List<Point<T>>>(size);
            for (int i = 0; i < size; i++)
                data_.Add(new List<Point<T>>());
        }



        public void Init(List<Point<T>> inputs)
        {
            lower_ = long.MaxValue;
            upper_ = 0;
            foreach (Point<T> it in inputs)
            {
                double val = it.Value();
                lower_ = Math.Min(lower_, val);
                upper_ = Math.Max(upper_, val);
            }
            lower_--;

            if (type_ == ToleranceBy.PPM)
                PPMInit(inputs);
            else if (type_ == ToleranceBy.Dalton)
                DaltonInit(inputs);
        }

        public bool Match(double expect, double baseValue)
        {
            int index = Index(expect);

            if (index < 0 || index >= (int)data_.Count)
                return false;

            if (data_[index].Count > 0)
                return true;

            int size = (int)data_.Count;
            if (index < size - 1)
            {
                var it = data_[index + 1];
                for (int i = 0; i < it.Count; i++)
                {
                    if (IsMatch(expect, it[i].Value(), baseValue))
                    {
                        return true;
                    }
                }
            }

            if (index > 0)
            {
                var it = data_[index - 1];
                int bin_size = it.Count;
                for (int i = bin_size - 1; i >= 0; i--)
                {
                    if (IsMatch(expect, it[i].Value(), baseValue))
                    {
                        return true;
                    }
                }
            }

            return false;
        }

        public bool Match(double expect)
        {
            return Match(expect, expect);
        }

        public List<Point<T>> Search(double expect)
        { return Search(expect, expect); }

        public List<Point<T>> Search(double expect, double baseValue)
        {
            List<Point<T>> result = new List<Point<T>>();
            int index = Index(expect);
            int size = (int)data_.Count;
            if (index < 0 || index >= size)
                return result;

            foreach (var it in data_[index])
            {
                result.Add(it);
            }


            if (index < size - 1)
            {
                var it = data_[index + 1];
                for (int i = 0; i < (int)it.Count; i++)
                {
                    if (IsMatch(expect, it[i].Value(), baseValue))
                    {
                        result.Add(it[i]);
                    }
                }
            }


            if (index > 0)
            {
                var it = data_[index - 1];
                int bin_size = (int)it.Count;
                for (int i = bin_size - 1; i >= 0; i--)
                {
                    if (IsMatch(expect, it[i].Value(), baseValue))
                    {
                        result.Add(it[i]);
                    }
                }
            }

            return result;
        }

        public List<T> SearchContent(double expect, double baseValue)
        {
            return Search(expect, baseValue).Select(r => r.Content()).ToList();
        }

        public List<T> SearchContent(double expect)
        {
            return SearchContent(expect, expect);
        }

        public void SetTolerance(double tol)
        {
            tolerance_ = tol;
        }

        public void SetToleranceBy(ToleranceBy by)
        {
            type_ = by;
        }

        public void Add(Point<T> point)
        {
            double expect = point.Value();
            int index = Index(expect);
            if (index >= 0 && index < (int)data_.Count)
                data_[index].Add(point);
        }

        bool IsMatch(double expect, double observe, double baseValue)
        {
            if (type_ == ToleranceBy.PPM)
            {
                return Math.Abs(expect - observe) / baseValue * 1000000.0 < tolerance_;
            }
            return Math.Abs(expect - observe) < tolerance_;
        }



        int Index(double expect)
        {
            if (type_ == ToleranceBy.Dalton)
                return (int)Math.Floor((expect - lower_) / tolerance_);
            double ratio = 1.0 / (1.0 - tolerance_ / 1000000);
            return (int)Math.Floor(Math.Log(expect * 1.0 / lower_) / Math.Log(ratio));
        }

        void DaltonInit(List<Point<T>> inputs)
        {
            DaltonInit();

            // assign the value
            foreach (Point<T> point in inputs)
            {
                double val = point.Value();
                if (val < lower_ || val > upper_)
                    continue;
                int index = (int)Math.Floor((val - lower_) / tolerance_);
                data_[index].Add(point);
            }
        }

        void PPMInit(List<Point<T>> inputs)
        {
            PPMInit();

            // assign the value
            double ratio = 1.0 / (1.0 - tolerance_ / 1000000);
            foreach (Point<T> point in inputs)
            {
                double val = point.Value();
                if (val < lower_ || val > upper_)
                    continue;
                int index = (int)Math.Floor(Math.Log(val / lower_) / Math.Log(ratio));
                data_[index].Add(point);
            }
        }

        public double Tolerance()
        {
            return tolerance_;
        }

        public ToleranceBy ToleranceType()
        {
            return type_;
        }
    }
}
