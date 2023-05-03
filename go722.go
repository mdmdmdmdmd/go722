package go722

import (
	"math"
)

type Options struct {
	Bitrate      Bitrate
	SampleRate8k bool
	Packed       bool
}

type Bitrate int

const (
	Bitrate48k Bitrate = iota + 6
	Bitrate56k
	Bitrate64k
)

type g722Band struct {
	s   int
	sp  int
	sz  int
	r   [3]int
	a   [3]int
	ap  [3]int
	p   [3]int
	d   [7]int
	b   [7]int
	bp  [7]int
	sg  [7]int
	nb  int
	det int
}

type context struct {
	ituTestMode   bool
	packed        bool
	eightK        bool
	bitsPerSample int
	x             [24]int
	band          [2]g722Band
	inBuffer      uint
	inBits        int
	outBuffer     uint
	outBits       int
}

var wl = [8]int{-60, -30, 58, 172, 334, 538, 1198, 3042}
var rl42 = [16]int{
	0, 7, 6, 5, 4, 3, 2, 1,
	7, 6, 5, 4, 3, 2, 1, 0,
}
var ilb = [32]int{
	2048, 2093, 2139, 2186, 2233, 2282, 2332, 2383,
	2435, 2489, 2543, 2599, 2656, 2714, 2774, 2834,
	2896, 2960, 3025, 3091, 3158, 3228, 3298, 3371,
	3444, 3520, 3597, 3676, 3756, 3838, 3922, 4008,
}
var qm4 = [16]int{
	0, -20456, -12896, -8968, -6288, -4240, -2584, -1200,
	20456, 12896, 8968, 6288, 4240, 2584, 1200, 0,
}
var qm2 = [4]int{-7408, -1616, 7408, 1616}
var qmfCoeffs = [12]int{
	3, -11, 12, 32, -210, 951, 3876, -805,
	362, -156, 53, -11,
}
var wh = [3]int{0, -214, 798}
var rh2 = [4]int{2, 1, 2, 1}

func ContextNew(options Options) (*context, error) {
	ctx := &context{}
	ctx.bitsPerSample = int(options.Bitrate)
	ctx.eightK = options.SampleRate8k
	if options.Packed && ctx.bitsPerSample != 8 {
		ctx.packed = true
	} else {
		ctx.packed = false
	}
	ctx.band[0].det = 32
	ctx.band[1].det = 8
	return ctx, nil
}

func saturate(amp int32) int16 {
	amp16 := int16(amp)
	if int32(amp16) == amp {
		return amp16
	}
	if amp > int32(math.MaxInt16) {
		return math.MaxInt16
	}
	return math.MinInt16
}

func (c *context) block4(index int, d int) {
	var wd1, wd2, wd3, i int

	// Block 4, RECONS
	c.band[index].d[0] = d
	c.band[index].r[0] = int(saturate(int32(c.band[index].s + d)))

	// Block 4, PARREC
	c.band[index].p[0] = int(saturate(int32(c.band[index].sz + d)))

	// Block 4, UPPOL2
	for i = 0; i < 3; i++ {
		c.band[index].sg[i] = c.band[index].p[i] >> 15
	}
	wd1 = int(saturate(int32(c.band[index].a[1] << 2)))

	if c.band[index].sg[0] == c.band[index].sg[1] {
		wd2 = -wd1
	} else {
		wd2 = wd1
	}

	if wd2 > 32767 {
		wd2 = 32767
	}

	if c.band[index].sg[0] == c.band[index].sg[2] {
		wd3 = (wd2 >> 7) + 128
	} else {
		wd3 = (wd2 >> 7) + -128
	}

	wd3 += (c.band[index].a[2] * 32512) >> 15

	if wd3 > 12288 {
		wd3 = 12288
	} else if wd3 < -12288 {
		wd3 = -12288
	}

	c.band[index].ap[2] = wd3

	// Block 4, UPPOL1
	c.band[index].sg[0] = c.band[index].p[0] >> 15
	c.band[index].sg[1] = c.band[index].p[1] >> 15
	if c.band[index].sg[0] == c.band[index].sg[1] {
		wd1 = 192
	} else {
		wd1 = -192
	}

	wd2 = (c.band[index].a[1] * 32640) >> 15

	c.band[index].ap[1] = int(saturate(int32(wd1 + wd2)))

	wd3 = int(saturate(int32(15360 - c.band[index].ap[2])))

	if c.band[index].ap[1] > wd3 {
		c.band[index].ap[1] = wd3
	} else if c.band[index].ap[1] < -wd3 {
		c.band[index].ap[1] = -wd3
	}

	// Block 4, UPZERO
	if d == 0 {
		wd1 = 0
	} else {
		wd1 = 128
	}

	c.band[index].sg[0] = d >> 15

	for i = 1; i < 7; i++ {
		c.band[index].sg[i] = c.band[index].d[i] >> 15

		if c.band[index].sg[i] == c.band[index].sg[0] {
			wd2 = wd1
		} else {
			wd2 = -wd1
		}

		wd3 = (c.band[index].b[i] * 32640) >> 15

		c.band[index].bp[i] = int(saturate(int32(wd2 + wd3)))
	}

	// Block 4, DELAYA
	for i = 6; i > 0; i-- {
		c.band[index].d[i] = c.band[index].d[i-1]
		c.band[index].b[i] = c.band[index].bp[i]
	}

	for i = 2; i > 0; i-- {
		c.band[index].r[i] = c.band[index].r[i-1]
		c.band[index].p[i] = c.band[index].p[i-1]
		c.band[index].a[i] = c.band[index].ap[i]
	}

	// Block 4, FILTEP
	wd1 = int(saturate(int32(c.band[index].r[1] + c.band[index].r[1])))
	wd1 = (c.band[index].a[1] * wd1) >> 15
	wd2 = int(saturate(int32(c.band[index].r[2] + c.band[index].r[2])))
	wd2 = (c.band[index].a[2] * wd2) >> 15
	c.band[index].sp = int(saturate(int32(wd1 + wd2)))

	// Block 4, FILTEZ
	c.band[index].sz = 0
	for i = 6; i > 0; i-- {
		wd1 = int(saturate(int32(c.band[index].d[i] + c.band[index].d[i])))
		c.band[index].sz += (c.band[index].b[i] * wd1) >> 15
	}
	c.band[index].sz = int(saturate(int32(c.band[index].sz)))

	// Block 4, PREDIC
	c.band[index].s = int(saturate(int32(c.band[index].sp + c.band[index].sz)))
}

func (c *context) Decode(data []uint8) []int16 {
	qm5 := [32]int{
		-280, -280, -23352, -17560, -14120, -11664, -9752, -8184,
		-6864, -5712, -4696, -3784, -2960, -2208, -1520, -880,
		23352, 17560, 14120, 11664, 9752, 8184, 6864, 5712,
		4696, 3784, 2960, 2208, 1520, 880, 280, -280,
	}
	qm6 := [64]int{
		-136, -136, -136, -136, -24808, -21904, -19008, -16704,
		-14984, -13512, -12280, -11192, -10232, -9360, -8576, -7856,
		-7192, -6576, -6000, -5456, -4944, -4464, -4008, -3576,
		-3168, -2776, -2400, -2032, -1688, -1360, -1040, -728,
		24808, 21904, 19008, 16704, 14984, 13512, 12280, 11192,
		10232, 9360, 8576, 7856, 7192, 6576, 6000, 5456,
		4944, 4464, 4008, 3576, 3168, 2776, 2400, 2032,
		1688, 1360, 1040, 728, 432, 136, -432, -136,
	}

	var dlowt, rlow, ihigh, dhigh, rhigh, xout1, xout2, wd1, wd2, wd3, code, outlen int
	dataLen := len(data)
	amp := make([]int16, dataLen*2)

	for j := 0; j < dataLen; {
		if c.packed {
			// Unpack the code bits
			if c.inBits < c.bitsPerSample {
				c.inBuffer |= uint(data[j]) << c.inBits
				j += 1
				c.inBits += 8
			}
			code = int(c.inBuffer & ((1 << c.bitsPerSample) - 1))
			c.inBuffer >>= c.bitsPerSample
			c.inBits -= c.bitsPerSample
		} else {
			code = int(data[j])
			j += 1
		}

		switch c.bitsPerSample {
		case 8:
			wd1 = code & 0x3F
			ihigh = (code >> 6) & 0x03
			wd2 = qm6[wd1]
			wd1 >>= 2
		case 7:
			wd1 = code & 0x1F
			ihigh = (code >> 5) & 0x03
			wd2 = qm5[wd1]
			wd1 >>= 1
		case 6:
			wd1 = code & 0x0F
			ihigh = (code >> 4) & 0x03
			wd2 = qm4[wd1]
		}

		// Block 5L, LOW BAND INVQBL
		wd2 = (c.band[0].det * wd2) >> 15
		// Block 5L, RECONS
		rlow = c.band[0].s + wd2
		// Block 6L, LIMIT
		if rlow > 16383 {
			rlow = 16383
		} else if rlow < -16384 {
			rlow = -16384
		}

		// Block 2L, INVQAL
		wd2 = qm4[wd1]
		dlowt = (c.band[0].det * wd2) >> 15

		// Block 3L, LOGSCL
		wd2 = rl42[wd1]
		wd1 = (c.band[0].nb * 127) >> 7
		wd1 += wl[wd2]
		if wd1 < 0 {
			wd1 = 0
		} else if wd1 > 18432 {
			wd1 = 18432
		}
		c.band[0].nb = wd1

		/* Block 3L, SCALEL */
		wd1 = (c.band[0].nb >> 6) & 31
		wd2 = 8 - (c.band[0].nb >> 11)
		if wd2 < 0 {
			wd3 = ilb[wd1] << -wd2
		} else {
			wd3 = ilb[wd1] >> wd2
		}
		c.band[0].det = wd3 << 2

		c.block4(0, dlowt)

		if !c.eightK {
			/* Block 2H, INVQAH */
			wd2 = qm2[ihigh]
			dhigh = (c.band[1].det * wd2) >> 15
			/* Block 5H, RECONS */
			rhigh = dhigh + c.band[1].s
			/* Block 6H, LIMIT */
			if rhigh > 16383 {
				rhigh = 16383
			} else if rhigh < -16384 {
				rhigh = -16384
			}

			/* Block 2H, INVQAH */
			wd2 = rh2[ihigh]
			wd1 = (c.band[1].nb * 127) >> 7
			wd1 += wh[wd2]
			if wd1 < 0 {
				wd1 = 0
			} else if wd1 > 22528 {
				wd1 = 22528
			}
			c.band[1].nb = wd1

			/* Block 3H, SCALEH */
			wd1 = (c.band[1].nb >> 6) & 31
			wd2 = 10 - (c.band[1].nb >> 11)
			if wd2 < 0 {
				wd3 = ilb[wd1] << -wd2
			} else {
				wd3 = ilb[wd1] >> wd2
			}
			c.band[1].det = wd3 << 2

			c.block4(1, dhigh)
		}

		if c.ituTestMode {
			amp[outlen] = int16(rlow << 1)
			outlen += 1
			amp[outlen] = int16(rhigh << 1)
			outlen += 1
		} else {
			if c.eightK {
				amp[outlen] = int16(rlow << 1)
				outlen += 1
			} else {
				// Apply the receive QMF
				for i := 0; i < 22; i++ {
					c.x[i] = c.x[i+2]
				}
				c.x[22] = rlow + rhigh
				c.x[23] = rlow - rhigh

				xout1 = 0
				xout2 = 0
				for i := 0; i < 12; i++ {
					xout2 += c.x[2*i] * qmfCoeffs[i]
					xout1 += c.x[2*i+1] * qmfCoeffs[11-i]
				}
				amp[outlen] = saturate(int32(xout1 >> 11))
				outlen += 1
				amp[outlen] = saturate(int32(xout2 >> 11))
				outlen += 1
			}
		}
	}
	return amp
}

func (c *context) Encode(data []int16) []uint8 {
	q6 := [32]int{
		0, 35, 72, 110, 150, 190, 233, 276,
		323, 370, 422, 473, 530, 587, 650, 714,
		786, 858, 940, 1023, 1121, 1219, 1339, 1458,
		1612, 1765, 1980, 2195, 2557, 2919, 0, 0,
	}
	iln := [32]int{
		0, 63, 62, 31, 30, 29, 28, 27,
		26, 25, 24, 23, 22, 21, 20, 19,
		18, 17, 16, 15, 14, 13, 12, 11,
		10, 9, 8, 7, 6, 5, 4, 0,
	}
	ilp := [32]int{
		0, 61, 60, 59, 58, 57, 56, 55,
		54, 53, 52, 51, 50, 49, 48, 47,
		46, 45, 44, 43, 42, 41, 40, 39,
		38, 37, 36, 35, 34, 33, 32, 0,
	}
	ihn := [3]int{0, 1, 0}
	ihp := [3]int{0, 3, 2}

	var dlow, dhigh, el, wd, wd1, ril, wd2, il4, ih2, wd3, eh, mih, i int
	/* Low and high band PCM from the QMF */
	var xlow, xhigh, g722Bytes int
	/* Even and odd tap accumulators */
	var sumeven, sumodd, ihigh, ilow, code int

	g722Bytes = 0
	xhigh = 0

	dataLen := len(data)
	output := make([]uint8, dataLen)

	for j := 0; j < dataLen; {
		if c.ituTestMode {
			xhigh = int(data[j] >> 1)
			xlow = xhigh
			j += 1
		} else {
			if c.eightK {
				xlow = int(data[j] >> 1)
				j += 1
			} else {
				/* Apply the transmit QMF */
				/* Shuffle the buffer down */
				for i := 0; i < 22; i++ {
					c.x[i] = c.x[i+2]
				}
				c.x[22] = int(data[j])
				j += 1
				c.x[23] = int(data[j])
				j += 1

				/* Discard every other QMF output */
				sumeven = 0
				sumodd = 0
				for i := 0; i < 12; i++ {
					sumodd += c.x[2*i] * qmfCoeffs[i]
					sumeven += c.x[2*i+1] * qmfCoeffs[11-i]
				}
				xlow = (sumeven + sumodd) >> 14
				xhigh = (sumeven - sumodd) >> 14
			}
		}
		/* Block 1L, SUBTRA */
		el = int(saturate(int32(xlow - c.band[0].s)))

		/* Block 1L, QUANTL */
		if el >= 0 {
			wd = el
		} else {
			wd = -(el + 1)
		}

		for i = 1; i < 30; i++ {
			wd1 = (q6[i] * c.band[0].det) >> 12
			if wd < wd1 {
				break
			}
		}
		if el < 0 {
			ilow = iln[i]
		} else {
			ilow = ilp[i]
		}

		/* Block 2L, INVQAL */
		ril = ilow >> 2
		wd2 = qm4[ril]
		dlow = (c.band[0].det * wd2) >> 15

		/* Block 3L, LOGSCL */
		il4 = rl42[ril]
		wd = (c.band[0].nb * 127) >> 7
		c.band[0].nb = wd + wl[il4]
		if c.band[0].nb < 0 {
			c.band[0].nb = 0
		} else if c.band[0].nb > 18432 {
			c.band[0].nb = 18432
		}

		/* Block 3L, SCALEL */
		wd1 = (c.band[0].nb >> 6) & 31
		wd2 = 8 - (c.band[0].nb >> 11)
		if wd2 < 0 {
			wd3 = (ilb[wd1] << -wd2)
		} else {
			wd3 = (ilb[wd1] >> wd2)
		}
		c.band[0].det = wd3 << 2

		c.block4(0, dlow)

		if c.eightK {
			/* Just leave the high bits as zero */
			code = (0xC0 | ilow) >> (8 - c.bitsPerSample)
		} else {
			/* Block 1H, SUBTRA */
			eh = int(saturate(int32(xhigh - c.band[1].s)))

			/* Block 1H, QUANTH */
			if eh >= 0 {
				wd = eh
			} else {
				wd = -(eh + 1)
			}
			wd1 = (564 * c.band[1].det) >> 12
			if wd >= wd1 {
				mih = 2
			} else {
				mih = 1
			}
			if eh < 0 {
				ihigh = ihn[mih]
			} else {
				ihigh = ihp[mih]
			}

			/* Block 2H, INVQAH */
			wd2 = qm2[ihigh]
			dhigh = (c.band[1].det * wd2) >> 15

			/* Block 3H, LOGSCH */
			ih2 = rh2[ihigh]
			wd = (c.band[1].nb * 127) >> 7
			c.band[1].nb = wd + wh[ih2]
			if c.band[1].nb < 0 {
				c.band[1].nb = 0
			} else if c.band[1].nb > 22528 {
				c.band[1].nb = 22528
			}

			/* Block 3H, SCALEH */
			wd1 = (c.band[1].nb >> 6) & 31
			wd2 = 10 - (c.band[1].nb >> 11)
			if wd2 < 0 {
				wd3 = (ilb[wd1] << -wd2)
			} else {
				wd3 = (ilb[wd1] >> wd2)
			}
			c.band[1].det = wd3 << 2

			c.block4(1, dhigh)
			code = ((ihigh << 6) | ilow) >> (8 - c.bitsPerSample)
		}

		if c.packed {
			/* Pack the code bits */
			c.outBuffer |= uint(code << c.outBits)
			c.outBits += c.bitsPerSample
			if c.outBits >= 8 {
				output[g722Bytes] = uint8(c.outBuffer & 0xFF)
				g722Bytes += 1
				c.outBits -= 8
				c.outBuffer >>= 8
			}
		} else {
			output[g722Bytes] = uint8(code)
			g722Bytes += 1
		}
	}
	return output
}
