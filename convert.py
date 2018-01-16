# -*- coding: utf-8 -*-
# KNMI HDF5 to Gamic HDF5 format converter
# Author - Kevin Broeren
# kevin@stormplatform.nl

import numpy as np
import os
import sys
import argparse
from datetime import datetime

try:
    import h5py
except:
    print "Error : Requires h5py"
    sys.exit()

parser = argparse.ArgumentParser(description='Provide HDF5 datafile.')
parser.add_argument('--knmi', help='HDF5 radar volume file')
args = parser.parse_args()
datafile = args.knmi


class Moment(object):
    name = ""
    unit = ""
    format = ""

    def __init__(self, name, unit, format):
        self.name = name
        self.unit = unit
        self.format = format


class CreateGamicHdf5Skeleton:
    """
    Creates a HDF5 skeleton in the Gamic HDF5 file format.
    """

    file = ""

    def __init__(self):
        self.createfile()

    def createfile(self):
        if os.path.isfile('tmp/temp.h5'):
            os.remove('tmp/temp.h5')
        self.file = h5py.File("tmp/temp.h5", "w")

    def createHowGroupWithAttributes(self, azimuth_beam, elevation_beam, sdp_manufacturer, site_name):
        how = self.file.create_group('how')
        how.attrs['azimuth_beam'] = azimuth_beam
        how.attrs['elevation_beam'] = elevation_beam
        how.attrs['sdp_manufacturer'] = sdp_manufacturer
        how.attrs['site_name'] = site_name

    def createWhatGroupWithAttributes(self, object, sets, sets_scheduled, version):
        what = self.file.create_group("what")
        what.attrs['object'] = object
        what.attrs['sets'] = sets
        what.attrs['sets_scheduled'] = sets_scheduled
        what.attrs['version'] = version

    def createWhereGroupWithAttributes(self, height, lat, lon):
        where = self.file.create_group('where')
        where.attrs['height'] = height
        where.attrs['lat'] = lat
        where.attrs['lon'] = lon

    def createScanHowGroupWithAttributes(self, scan, timestamp, bin_count, range_start, scan_speed, range, range_step,
                                         range_samples, PRF, angle_sync, time_samples, angle_step, unfolding,long_pulse,
                                         clutter_filter_number, radar_wave_length, output, output64, half_resolution,
                                         ray_count, azi_start, azi_stop, ele_start, ele_stop, ele_step, azimuth,
                                         elevation, pulse_width_us, malfunc, raster):
        how = scan.create_group('how')
        how.attrs['timestamp'] = timestamp
        how.attrs['bin_count'] = bin_count
        how.attrs['range_start'] = range_start
        how.attrs['scan_speed'] = scan_speed
        how.attrs['range'] = range
        how.attrs['range_step'] = range_step
        how.attrs['range_samples'] = range_samples
        how.attrs['PRF'] = PRF
        how.attrs['angle_sync'] = angle_sync
        how.attrs['time_samples'] = time_samples
        how.attrs['angle_step'] = angle_step
        how.attrs['unfolding'] = unfolding
        how.attrs['long_pulse'] = long_pulse
        how.attrs['clutter_filter_number'] = clutter_filter_number
        how.attrs['radar_wave_length'] = radar_wave_length
        how.attrs['output'] = output
        how.attrs['output64'] = output64
        how.attrs['half_resolution'] = half_resolution
        how.attrs['ray_count'] = ray_count
        how.attrs['azi_start'] = azi_start
        how.attrs['azi_stop'] = azi_stop
        how.attrs['ele_start'] = ele_start
        how.attrs['ele_stop'] = ele_stop
        how.attrs['ele_step'] = ele_step
        how.attrs['azimuth'] = azimuth
        how.attrs['elevation'] = elevation
        how.attrs['pulse_width_us'] = pulse_width_us
        how.attrs['malfunc'] = malfunc
        how.attrs['raster'] = raster
        bite = how.create_group('bite')
        bite.create_group('raw')
        how.create_group('extended')
        how.create_group('raw_sdp_param')

    def createScanWhatGroupWithAttributes(self, scan, product, scan_type, descriptor_count):
        what = scan.create_group('what')
        what.attrs['product'] = product
        what.attrs['scan_type'] = scan_type
        what.attrs['descriptor_count'] = descriptor_count

    def getFile(self):
        return self.file


class KNMI_Converter:
    """
    Converts the contents of an KNMI HDF5 file into Gamic HDF5 file format. 
    """

    def __init__(self, input_file):

        volumefile = h5py.File(input_file, "r")

        skeleton = CreateGamicHdf5Skeleton()

        skeleton.createHowGroupWithAttributes(
            azimuth_beam=1.0,
            elevation_beam=1.0,
            sdp_manufacturer='Stormplatform',
            site_name=volumefile['radar1'].attrs['radar_location']
        )

        skeleton.createWhatGroupWithAttributes(
            object='PVOL',
            sets=volumefile['overview'].attrs['number_scan_groups'],
            sets_scheduled=volumefile['overview'].attrs['number_scan_groups'],
            version=9
        )

        skeleton.createWhereGroupWithAttributes(
            height=55,
            lat=volumefile['radar1'].attrs['radar_location'][1],
            lon=volumefile['radar1'].attrs['radar_location'][0]
        )

        newhdf5file = skeleton.getFile()

        # Loop the scans
        number_scan_groups = volumefile['overview'].attrs['number_scan_groups']
        for scan_number in range(0, number_scan_groups):
            scan_group_name = 'scan' + str(scan_number)
            knmi_scan_group_name = 'scan' + str(scan_number + 1)

            # create Scan group
            scangroup = newhdf5file.create_group(scan_group_name)
            skeleton.createScanHowGroupWithAttributes(
                scan=scangroup,
                timestamp=self.convertKnmiTimestamp(volumefile[knmi_scan_group_name].attrs['scan_datetime']),
                bin_count='',
                range_start=0.0,
                scan_speed=volumefile[knmi_scan_group_name].attrs['scan_antenna_velocity'],
                range=volumefile[knmi_scan_group_name].attrs['scan_number_range']+000,  # In meters
                range_step='',
                range_samples=1,
                PRF=volumefile[knmi_scan_group_name].attrs['scan_high_PRF'],
                angle_sync='',
                time_samples='',
                angle_step='',
                unfolding='',
                long_pulse='',
                clutter_filter_number='',
                radar_wave_length='',
                output='',
                output64='',
                half_resolution='',
                ray_count='',
                azi_start=volumefile[knmi_scan_group_name].attrs['scan_start_azim'],
                azi_stop='',
                ele_start=volumefile[knmi_scan_group_name].attrs['scan_number_range'],
                ele_stop='',
                ele_step='',
                azimuth=volumefile[knmi_scan_group_name].attrs['scan_number_azim'],
                elevation=volumefile[knmi_scan_group_name].attrs['scan_elevation'],
                pulse_width_us='',
                malfunc='',
                raster=''
            )

            skeleton.createScanWhatGroupWithAttributes(
                scan=scangroup,
                product='SCAN',
                scan_type='PPI',
                descriptor_count=18
            )

            # Loop the moments in the scan groups
            for idx, scantype in enumerate(volumefile[knmi_scan_group_name]):
                moment_group_name = 'moment_' + str(idx - 1)
                moment_data = volumefile[knmi_scan_group_name][scantype]
                moment_data_array = np.array(moment_data)

                # Only get the 2 dimensional data arrays
                if moment_data_array.ndim == 2:
                    # Create new dataset in the scan group
                    scan_moment_dataset = scangroup.create_dataset(moment_group_name,
                                                                   (moment_data_array.shape[0],
                                                                    moment_data_array.shape[1]),
                                                                   chunks=(moment_data_array.shape[0],
                                                                           moment_data_array.shape[1]),
                                                                   compression="gzip",
                                                                   compression_opts=6,
                                                                   data=moment_data_array)

                    # Get the moment attributes and add to the dataset
                    moment_attributes = self.KnmiScanTypeToMoment(scantype)
                    for attrname, attrvalue in moment_attributes.__dict__.iteritems():
                        scan_moment_dataset.attrs[attrname] = attrvalue

    def KnmiScanTypeToMoment(self, scantype):
        switcher = {
            'scan_CCOR_data': Moment(name='CCOR', unit='dB', format='UV16'),
            'scan_CCORv_data': Moment(name='CCORv', unit='dB', format='UV16'),
            'scan_KDP_data': Moment(name='KDP', unit='°/km', format='UV16'),
            'scan_PhiDP_data': Moment(name='PHIDP', unit='°', format='UV16'),
            'scan_RhoHV_data': Moment(name='RHOHV', unit='Coefficient', format='UV16'),
            'scan_SQI_data': Moment(name='SQI', unit='Coefficient', format='UV16'),
            'scan_SQIv_data': Moment(name='SQIv', unit='Coefficient', format='UV16'),
            'scan_V_data': Moment(name='V', unit='m/s', format='UV16'),
            'scan_Vv_data': Moment(name='Vv', unit='m/s', format='UV16'),
            'scan_W_data': Moment(name='W', unit='m/s', format='UV16'),
            'scan_Wv_data': Moment(name='Wv', unit='m/s', format='UV16'),
            'scan_Z_data': Moment(name='Z', unit='dBZ', format='UV16'),
            'scan_Zv_data': Moment(name='Zv', unit='dBZ', format='UV16'),
            'scan_uPhiDP_data': Moment(name='UPHIDP', unit='°', format='UV16'),
            'scan_uZ_data': Moment(name='UZ', unit='dBZ', format='UV16'),
            'scan_uZv_data': Moment(name='UZv', unit='dBZ', format='UV16'),
            'scan_CPA_data': Moment(name='CPA', unit='na', format='UV16'),
            'scan_CPAv_data': Moment(name='CPAv', unit='na', format='UV16'),
            'scan_TX_power': Moment(name='TX', unit='na', format='UV32'),
        }
        return switcher.get(scantype, "invalid scan type")

    def convertKnmiTimestamp(self, timestamp):
        # input template: 22-JUN-2017;14:25:04.000
        # output template: 2017-07-22T14:25:04Z
        d = datetime.strptime(timestamp, '%d-%b-%Y;%H:%M:%S.000')
        return d.strftime('%Y-%m-%dT%H:%M:%SZ')

# Convert the file
KNMI_Converter(datafile)
