from datetime import datetime

class SolarLookupTime(object):
    x = [-10, -15, -20, -25, -30, -35, -40, -44]
    allhours = [18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
    def __init__(self,name,ymd1,ymd2,obstime,hours,latitude=x):
        self.name = name
        self.d1 = datetime.strptime(str(ymd1), "%Y%m%d")
        self.d2 = datetime.strptime(str(ymd2), "%Y%m%d")
        self.obstime = obstime
        self.hours = hours
        if self.hours is None: self.hours = self.allhours
        self.lats = latitude
    def __repr__(self):
        d = {}
        d['name'] = self.name
        d['start_date'] = self.d1.strftime("%Y-%m-%d")
        d['end_date'] = self.d2.strftime("%Y-%m-%d")
        d['latitudes'] = self.lats
        d['observation_time'] = self.obstime
        d['applicable_hours'] = self.hours
        return str(d)
    def __str__(self):
        return self.__repr__()
    def is_applicable(self,datehour):
        if datehour.date() >= self.d1.date() and \
           datehour.date() <= self.d2.date() and \
           datehour.hour in self.hours: return True
        return False
    def get_lookup(self):
        return self.lats, self.obstime


def get_lookup_list():
    """
    The lookup table defines the minutes after a start hour at
    five degree increments of latitude from -10 to -44 N (defined as x in
    the SolarLookupTime class).
    """
    # Choose lookup table based on date and hour
    lookup = []
    hours_a = [18, 19, 20, 21, 23, 0, 1, 2, 3, 5, 6, 7, 8, 9, 11]
    hours_b = [22, 4, 10]
    lookup.append(SolarLookupTime('GMS-4 A',19891231,19921231,
                  [45.7, 46.7, 47.7, 48.7, 49.6, 50.5, 51.2, 51.8],hours_a))
    lookup.append(SolarLookupTime('GMS-4 B',19891231,19921231,
                  [38.7, 39.7, 40.7, 41.7, 42.6, 43.5, 44.2, 44.8],hours_b))
    lookup.append(SolarLookupTime('GMS-4 A',19930101,19940630,
                  [47.2, 48.2, 49.3, 50.2, 51.1, 52, 52.7, 53.3],hours_a))
    lookup.append(SolarLookupTime('GMS-4 B',19930101,19940630,
                  [40.7, 41.7, 42.8, 43.7, 44.6, 45.5, 46.2, 46.8],hours_b))
    lookup.append(SolarLookupTime('GMS-4 A',19940701,19950610,
                  [46.7, 47.7, 48.8, 49.7, 50.6, 51.5, 52.2, 52.8],hours_a))
    lookup.append(SolarLookupTime('GMS-4 B',19940701,19950610,
                  [40.5, 41.5, 42.6, 43.5, 44.4, 45.3, 46, 46.6],hours_b))
    lookup.append(SolarLookupTime('GMS-5 A',19950611,20030520,
                  [46.7, 47.7, 48.8, 49.7, 50.6, 51.5, 52.2, 52.8],hours_a))
    lookup.append(SolarLookupTime('GMS-5 B',19950611,20030520,
                  [39.7, 40.7, 41.8, 42.7, 43.6, 44.5, 45.2, 45.8],hours_b))
    lookup.append(SolarLookupTime('GEOS-9 A',20030521,20051031,
                  [39.9, 41, 42, 43, 43.9, 44.7, 45.5, 46],hours_a))
    lookup.append(SolarLookupTime('GEOS-9 B',20030521,20051031,
                  [27.9, 29, 30, 31, 31.9, 32.7, 33.5, 34],hours_b))
    lookup.append(SolarLookupTime('MTSAT-1R',20051101,20100930,
                  [46.2, 47.2, 48.3, 49.2, 50.1, 51, 51.7, 52.3],None))
    lookup.append(SolarLookupTime('MTSAT-2',20100701,20120930,
                  [44.7, 45.7, 46.8, 47.7, 48.6, 49.5, 50.2, 50.8],None))
    return lookup
