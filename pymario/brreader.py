import struct
import os

__all__ = ['BRReader']

# TODO: handle ns5


def to_pystr(sz):
    """convert NULL terminated string to Python string"""
    return sz.split('\0')[0]


class BufferedFile:
    def __init__(self, fn, mode):
        assert mode == 'rb'
        self.fp = None
        self.close()
        self.fp = open(fn, 'rb')

    def close(self):
        if self.fp != None: self.fp.close()
        self.pos = 0
        self.fp_pos = 0
        self.buf = ''
        self.buf_off = -1
        self.fp = None

    def seek(self, offset, whence=0):
        if whence == 0: self.pos = offset
        elif whence == 1: self.pos += offset
        else: raise ValueError, 'Not implemented'

    def tell(self):
        return self.pos

    def read(self, size, readahead=None):
        if readahead is None:
            return self._read(size)

        # if data fetch is needed..
        if (self.buf_off > self.pos) or (self.buf_off + len(self.buf) < self.pos + size):
            self.buf_off = self.pos
            self.buf = self._read(readahead)
            self.pos = self.buf_off

        ib = self.pos - self.buf_off
        ie = ib + size
        r = self.buf[ib:ie]
        self.pos += len(r)
        return r

    def _read(self, size):
        f = self.fp
        # to reduce unnecessary fseek
        if self.fp_pos != self.pos:
            self.fp_pos = self.pos
            f.seek(self.pos)
        r = f.read(size)
        new_pos = self.pos + len(r)
        self.pos = new_pos
        self.fp_pos = new_pos
        return r


# BlackRockReader class
class BRReader:
    # Properties -----------------------------------------------------------
    MAGIC_NEV = 'NEURALEV'
    MAGIC_NSX = 'NEURALCD'   # .nsx file magic header
    MAGIC_NSX_EXT = 'CC'     # extended block magic header in .nsx
    L_EXTBLOCK = 32
    L_EXTID = 8
    L_SIGHEADER = 8
    L_NSX_BCHAN = 2          # number of bytes/channel/data pt (currently 2)
    I_NSX_BASEFREQ = 30000   # base frequency of maximum sampling rate in Hz

    @property
    def valid(self):
        return self._valid
    @property
    def t_sample(self):
        return self._t_sample
    @property
    def chn_info(self):
        return self._chn_info.copy()
    @property
    def spike_id(self):
        # returns list of electrode id that has spiking info
        return self._spike_id[:]
    @property
    def nsx_chn_order(self):
        return self._nsx_chn_order[:]

    # Methods -------------------------------------------------------------
    def __init__(self, nev_filename=None, nsx_filename=None):
        self.set_files(nev_filename=nev_filename, nsx_filename=nsx_filename)
        self._nev_fp = None
        self._nsx_fp = None
        self._reset()


    def _reset(self):
        if self._nev_fp != None:
            self._nev_fp.close()
        if self._nsx_fp != None:
            self._nsx_fp.close()

        self._nev_fp = None
        self._nsx_fp = None
        self._valid = False
        self._i_vmajor = -1
        self._i_vminor = -1
        self._i_nev_flags = -1
        self._l_header = -1
        self._l_packet = -1
        self._t_res = -1
        self._t_sample = -1
        self._n_extheader = -1
        self._chn_info = {}
        self._s_origin = None
        self._s_app = ''
        self._s_comment = ''
        self._p_databegin = -1
        self._spike_id = None

        self._nsx_l_header = -1
        self._nsx_s_label = ''
        self._nsx_comment = ''
        self._nsx_t_res = -1
        self._nsx_t_period = -1
        self._nsx_n_chan = -1
        self._nsx_chn_info = None
        self._nsx_valid = False
        self._nsx_index = None
        self._nsx_chn_order = None
        self._nsx_i_res = -1
        self._nsx_i_period = -1


    def set_files(self, nev_filename=None, nsx_filename=None):
        self.nev_filename = nev_filename
        self.nsx_filename = nsx_filename


    def close(self):
        self._reset()


    def open(self):
        return self._open_nev() and self._open_nsx()

    def _open_nev(self):
        """.nev file header reader"""
        nev_filename = self.nev_filename
        if nev_filename is None: return True
        #f = open(nev_filename, 'rb')
        f = BufferedFile(nev_filename, 'rb')

        # header and packet size
        s_header, i_vmajor, i_vminor, i_flags = struct.unpack('<8sBBH', f.read(12))
        if s_header != BRReader.MAGIC_NEV:
            raise IOError, 'Bad magic bytes'
        if i_vmajor != 2 or i_vminor != 2:
            raise IOError, 'Incompatible version (only works with v2.2)'
        l_header, l_packet = struct.unpack('<LL', f.read(8))
        l_extblock = BRReader.L_EXTBLOCK
        l_extid = BRReader.L_EXTID
        # t_res: time resolution of time stamp.
        #       (the frequency (counts per second) of the global clock used to index
        #       the time samples of the individual data packet entries.)
        # t_sample: time resolution of wav sampling
        t_res, t_sample = struct.unpack('<LL', f.read(8))    # in times per second (Hz)
        t_res = 1000000. / float(t_res)         # converting to us
        t_sample = 1000000. / float(t_sample)   # converting to us
        s_t_origin = f.read(16)          # time origin in Win32 SYSTEMTIME struct
        s_app = to_pystr(f.read(32))     # application info
        s_comment = to_pystr(f.read(256))
        # num of ext headers
        n_extheader = struct.unpack('<L', f.read(4))[0]
        chn_info = {}
        spike_id = []

        # read all ext headers --
        for i in range(n_extheader):
            block = f.read(l_extblock)
            blockid = block[:l_extid]

            if blockid == 'NEUEVWAV':
                elec_id, phy_conn, conn_pin, dig_factor, energy_thr, high_thr, \
                        low_thr, n_sorted, n_bwav = \
                            struct.unpack_from('<HBBHHhhBB', block[l_extid:])
                # dig_nv: digitization factor in nV
                chn_info[elec_id] = {'dig_nv':dig_factor, 'n_bwav':n_bwav,
                                     'energy_thr':energy_thr,
                                     'high_thr':high_thr, 'low_thr':low_thr} 
                spike_id.append(elec_id)

        # done.
        statinfo = os.stat(nev_filename)
        self._nev_fp = f
        self._i_vmajor = i_vmajor
        self._i_vminor = i_vminor
        self._i_nev_flags = i_flags
        self._valid = True
        self._l_header = l_header
        self._l_packet = l_packet
        self._n_packets = int(round(float(statinfo.st_size - l_header) / l_packet))
        self._t_res = t_res
        self._t_sample = t_sample
        self._n_extheader = n_extheader
        self._chn_info = chn_info
        self._s_origin = s_t_origin
        self._s_app = s_app
        self._s_comment = s_comment
        self._p_databegin = f.tell()
        self._spike_id = spike_id

        return True


    def _open_nsx(self):
        """.nsx file header reader"""
        nsx_filename = self.nsx_filename
        if nsx_filename is None: return True
        f = open(nsx_filename, 'rb')

        # Magic strig, major version, minor version, #bytes of all header info
        s_magic, i_vmajor, i_vminor, l_header = struct.unpack('<8sBBI', f.read(14))
        if s_magic != BRReader.MAGIC_NSX:
            raise ValueError, 'Bad NSX magic bytes'
        if i_vmajor != 2 or i_vminor != 2:
            raise ValueError, 'Incompatible version (only works with v2.2)'
        
        # label of the sampling group, comment 
        s_label = to_pystr(f.read(16))
        s_comment = to_pystr(f.read(256))

        # i_period: Number of 1/30,000 seconds between data points 
        #       (e.g., sampling rate of 30 kS/s = 1; 10 kS/s = 3)
        # i_res: frequency (counts per second) of the global clock used to index
        #       the time samples
        i_period, i_res = struct.unpack('<LL', f.read(8)) 
        t_res = 1000000. / float(i_res)                # converting to us
        t_period = i_period * 1. / BRReader.I_NSX_BASEFREQ * 1000000.   # converting to us
        s_t_origin = f.read(16)          # time origin in Win32 SYSTEMTIME struct
        # n_chan: number of channels per data point. (should be equal to the
        #         number of extended headers)
        n_chan = struct.unpack('<L', f.read(4))[0]
        
        chn_info = {}
        chn_order = []   # the actual data will be in this order
        # read all extended headers
        for i in range(n_chan):
            s_magic = f.read(2)
            if s_magic != BRReader.MAGIC_NSX_EXT:
                raise ValueError, 'Bad NSX magic bytes in extended header'

            # electrode id, electrode label, physical connector (e.g., bank A is 1)
            #    connector pin number in the bank
            i_elecid, s_eleclbl, i_conn, i_pin = \
                    struct.unpack('<H16sBB', f.read(20))
            s_eleclbl = to_pystr(s_eleclbl)
            
            # min/max of digital value
            i_mindig, i_maxdig = struct.unpack('<hh', f.read(4))
            # min/max of analog value. unit is determined by the next field.
            i_minalg, i_maxalg = struct.unpack('<hh', f.read(4))
            # actual units for i_minalg, i_maxalg. in string.
            s_algunit = to_pystr(f.read(16))

            # i_highf_cutoff: highpass filter freq cutoff in mHz
            # i_highf_order: Order of the filter used for high frequency cutoff:
            #   0 = NONE
            # i_highf_type: type of filter used for high frequency cutoff:
            #   0 = NONE, 1 = Butterworth
            i_highf_cutoff, i_highf_order, i_highf_type = struct.unpack('<LLH', f.read(10))
            # ... and everything is same for low freq cutoff
            i_lowf_cutoff, i_lowf_order, i_lowf_type = struct.unpack('<LLH', f.read(10))

            chn_info[i_elecid] = {'eleclbl': s_eleclbl, 'connector': i_conn,
                                  'pin': i_pin, 'mindig': i_mindig, 'maxdig':
                                  i_maxdig, 'minalg': i_minalg, 'maxalg':
                                  i_maxalg, 'algunit': s_algunit,
                                  'highf_cutoff': i_highf_cutoff,
                                  'highf_order': i_highf_order, 'highf_type': i_highf_type,
                                  'lowf_cutoff': i_lowf_cutoff,
                                  'lowf_order': i_lowf_order, 'lowf_type': i_lowf_type}
            chn_order.append(i_elecid)

        # nsx data index --------------------------------
        nsx_index = []
        n_datpt_cumul = 0
        while True:
            s_header = f.read(1)
            # header of a data packet should be 0x01
            if s_header != '\x01':
                if s_header == '': 
                    break    # reached EOF. finishing...
                else:
                    raise ValueError, 'Bad data packet header (got: %s)' % s_header

            # time stamp at the beginning, # data points in this packet
            t_timestamp, n_datpt = struct.unpack('<LL', f.read(8))
            # number of bytes/channel/data pt (currently 2)
            n_bchan = BRReader.L_NSX_BCHAN
            # start position of the data points
            pos = f.tell()
            # save the entry
            nsx_index.append({'ts_begin': t_timestamp * t_res,                              # onset of the first data point in this packet
                              'ts_end': (t_timestamp + (n_datpt - 1) * i_period) * t_res,   # onset of the last data point in this packet
                              'ts_endblk': (t_timestamp + (n_datpt) * i_period) * t_res,    # timestamp of the end of this packet
                              'n_datpt': n_datpt,
                              'pt_begin': n_datpt_cumul,
                              'pt_end': n_datpt_cumul + (n_datpt - 1),
                              'pos': pos})
            # move to the next entry
            f.seek(n_bchan * n_datpt * n_chan, 1)
            n_datpt_cumul += n_datpt

        self._nsx_fp = f
        self._nsx_l_header = l_header
        self._nsx_s_label = s_label
        self._nsx_comment = s_comment
        self._nsx_t_res = t_res
        self._nsx_i_res = i_res
        self._nsx_t_period = t_period
        self._nsx_i_period = i_period
        self._nsx_n_chan = n_chan
        self._nsx_chn_info = chn_info
        self._nsx_valid = True
        self._nsx_index = nsx_index
        self._nsx_chn_order = chn_order

        # move to the first data point
        self.nsx_goto_first_data()
        # done.
        return True


    # .nev file methods ====================================================================
    def read_once(self, pos=None, proc_wav=False, readahead=None, readahead_abs=None):
        """read one particular block: must be on the right position!!!"""
        l_sigheader = BRReader.L_SIGHEADER   # signal block header size
        l_packet = self._l_packet            # signal block size
        chn_info = self._chn_info
        t_res = self._t_res
        f = self._nev_fp

        if readahead != None: readahead_abs = (readahead + 1) * l_packet
        if pos != None: f.seek(pos)
        else: pos = f.tell()
        signal_block = f.read(l_packet, readahead=readahead_abs)
        if signal_block == '': return None   # reached EOF

        sig_t0, sig_id, ext_info = struct.unpack('<LHBx', signal_block[:l_sigheader])
        sig_t = sig_t0 * t_res   # converting to us
        signal_body = signal_block[l_sigheader:]
        all_data = {'id':sig_id, 'timestamp':sig_t, 'file_pos':pos}

        # if it is digital input information:
        if sig_id == 0:
            # TODO: handle other info (e.g., serial, etc.)
            dat_low, dat_high = struct.unpack_from('BB', signal_body)
            all_data['data_low8'] = dat_low
            all_data['data'] = dat_high

        # otherwise, it's waveform:
        else:             
            # for now, if proc_wav == False, nothing is done.
            elec_id = sig_id
            if proc_wav:
                n_bwav = chn_info[elec_id]['n_bwav']
                n_point = (l_packet - l_sigheader) / n_bwav

                if n_bwav == 1: fmt = '<' + 'b'*n_point
                elif n_bwav == 2: fmt = '<' + 'h'*n_point
                else: print 'Error!'

                wav = struct.unpack(fmt, signal_body)   # data points (in uV)
                all_data['waveform'] = wav

        return all_data


    def goto_first_data(self):
        self._nev_fp.seek(self._p_databegin)   # to the beginning of data

    def extract_all(self, proc_wav=False, only_timestamp=False):
        """one-shot nev file reader. do not use for big files."""
        if not self._valid:
            return False
        f = self._nev_fp
        self.goto_first_data()

        # read all data signal blocks ----------
        t_timestamp = []    # when the NSP got the time stamps?  (in us)
        d_timestamp = []    # what is the actual value of the stamps?
        i_timestamp = []    # index in "whole_signal"
        all_signal = []     # [(timestamp, data), ...] 
        l_sigheader = BRReader.L_SIGHEADER   # signal block header size

        while True:
            the_data = self.read_once(proc_wav=proc_wav)
            if the_data == None: break       # reached EOF

            # if it is digital input information:
            if the_data['id'] == 0:
                # TODO: handle other info (e.g., serial, etc.)
                sig_t, dat_low, dat_high = \
                        the_data['timestamp'], the_data['data_low8'], the_data['data']

                t_timestamp.append(sig_t)
                d_timestamp.append(dat_high)
                all_signal.append(the_data)
                i_timestamp.append(len(all_signal) - 1)

            # otherwise, it's waveform:
            elif not only_timestamp:             
                all_signal.append(the_data)

        return t_timestamp, d_timestamp, i_timestamp, all_signal


    # .nsx file methods =============================================================================================
    def _nsx_get_pt_offset(self, i_pt, full=False):
        """get the `i_pt`-th (`i_pt`: 0-based) data point's offest (and the offset of the last pt in this packet, etc.)"""
        idx_remainder = i_pt
        # number of bytes/channel/data pt (currently 2)
        n_bchan = BRReader.L_NSX_BCHAN
        # number of channels
        n_chan = self._nsx_n_chan
        found = False

        for i_packet, idx_info in enumerate(self._nsx_index):
            n_datpt = idx_info['n_datpt']
            if idx_remainder < n_datpt:
                found = True
                # pos: file offset for the `i_pt`-th data point
                # end: file offset for the last data point in this packet
                pos = idx_info['pos'] + n_bchan * n_chan * idx_remainder
                end = idx_info['pos'] + n_bchan * n_chan * (n_datpt - 1)
                n_remain = n_datpt - idx_remainder 
                break
            idx_remainder -= n_datpt

        if not found:
            raise IOError, 'Error: out of bound'

        # if we need just the file offset:
        if not full: return pos

        # if we want full information:
        # end: the file offset of the last pt in this packet
        # n_remain: number of remaining data points in this packet (including the `i_pt`-th data point)
        # i_next: the first data point # of the next block
        if i_packet == len(self._nsx_index) - 1: i_next = None
        else: i_next = self._nsx_index['pt_begin']
        return pos, end, n_remain, i_next


    def _nsx_get_ts_pt(self, ts):
        """get the nearest data point # for the specified timestamp"""
        # global timestamp resolution
        t_res = self._nsx_t_res
        i_period = self._nsx_i_period
        nsx_index = self._nsx_index
        ts = float(ts)
        found = False

        for idx_info in nsx_index:
            if idx_info['ts_begin'] <= ts and ts <= idx_info['ts_end']:
                found = True
                break
        if not found:
            ts_begins = [-float('inf')] + [info['ts_endblk'] for info in nsx_index]
            ts_ends = [info['ts_begin'] for info in nsx_index] + [float('inf')]
            for i, (begin, end) in enumerate(zip(ts_begins, ts_ends)):
                if not (begin <= ts and ts <= end): continue

                if end - ts < ts - begin:   # if near to the ts_end,
                    return nsx_index[i]['pt_begin']
                else:   # otherwise
                    return nsx_index[i - 1]['pt_end']
            # still not found. why? god whyyyy?!
            raise IOError, 'Error: out of bound - which is almost impossible (bug?)'

        # ts is somewhere in the packet described by `idx_info`
        t_base = idx_info['ts_begin']
        n_rel = int(round((ts - t_base) / t_res / i_period))
        i_pt = max(idx_info['pt_begin'], idx_info['pt_begin'] + n_rel)
        i_pt = min(idx_info['pt_end'], i_pt)
        return i_pt


    # some .nsx readers ------------------------------------------
    def nsx_read_once(self, pos=None, data_only=True, list_form=False):
        """read the current data point (for all channels). note: there's no safety check mechanism!!!"""
        f = self._nsx_fp
        chn_order = self._nsx_chn_order

        if pos != None: f.seek(pos)
        else: pos = f.tell()

        # number of bytes/channel/data pt (currently 2)
        n_bchan = BRReader.L_NSX_BCHAN
        # number of channels
        n_chan = self._nsx_n_chan
        raw = f.read(n_bchan * n_chan)
        if raw == '': return None
        fmt = '<' + 'h' * n_chan
        data = struct.unpack(fmt, raw)
        
        if list_form:
            datapt = list(data)
        else:
            datapt = {}
            for ch, pt in zip(chn_order, data): datapt[ch] = pt

        # if we just want data points only
        if data_only: return datapt
        # otherwise...
        return {'file_pos':pos, 'datapt':datapt}


    def nsx_goto_first_data(self):
        """rewind to the first data point"""
        pos = self._nsx_get_pt_offset(0)
        self._nsx_fp.seek(pos)

    def nsx_read_one_pt(self, i_pt, data_only=True, list_form=False):
        """read the `i_pt`-th (`i_pt`: 0-based) data point"""
        if not self._nsx_valid:
            raise IOError, '.nsx file not ready'

        pos = self._nsx_get_pt_offset(i_pt)
        return self.nsx_read_once(pos=pos, data_only=data_only,
                                  list_form=list_form)


    def nsx_read_range_ts(self, ts_begin, ts_end, 
                          data_only=True, list_form=False):
        if ts_end < ts_begin: raise ValueError, 'Bad range: ts_begin > ts_end'
        i_pt_begin = self._nsx_get_ts_pt(ts_begin)
        i_pt_end = self._nsx_get_ts_pt(ts_end)
        return self.nsx_read_range_pt(i_pt_begin, i_pt_end, data_only=data_only,
                                     list_form=list_form)


    def nsx_read_range_pt(self, i_pt_begin, i_pt_end, 
                          data_only=True, list_form=False):
        if i_pt_end < i_pt_begin: raise ValueError, 'Bad range: i_pt_begin > i_pt_end'
        # read the data points in the [i_pt_begin, i_pt_end] packet by packet
        while True:
            # read the information of this packet
            pos, end, n_remain, i_next = self._nsx_get_pt_offset(i_pt_begin, full=True)
            self._nsx_fp.seek(pos)
            i_pt = i_pt_begin
            for _ in range(min(n_remain, i_pt_end - i_pt_begin + 1)):
                yield self.nsx_read_once(data_only=data_only,
                                         list_form=list_form)
                i_pt =+ 1
            if i_next == None or i_pt > i_pt_end: break
            i_pt_begin = i_next
