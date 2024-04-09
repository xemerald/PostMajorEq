# PostMajorEq
The Major Earthquake waveform record processing program for SAC, TANK or miniSEED files.

## Dependencies
Nothing special, it can be compiled under only standard C library.

## Supported Platforms
- Linux
- MacOS (**Not pass testing yet**)

## Building & Installation
- Linux/MacOS
	- Simply run `make`
	- Then you can run `sudo make install` to install this program to `/usr/local/bin`

## Usage
- `postmajor -h` show some helping tips.
- `postmajor -v` show version information.
- `postmajor -f SAC <input eq. info> <input station list> <input seismic data>` or `postmajor <input eq. info> <input station list> <input seismic data>` process the input **SAC** format file(s) & output the result to the standard output.
- `postmajor -f MSEED <input eq. info> <input station list> <input seismic data>` process the input **miniSEED** format file(s).
- `postmajor <input eq. info> <input station list> <input seismic data> > <output path>` process the input **SAC** format file(s) & redirect the result to the output path.
- `postmajor -c <input eq. info> <input station list> <input seismic data>` process the input **SAC** format file(s) & append the station coordinate to the result.

## Earthquake information & Station list file content
Please refer to the example files.

## Output field description
- Default:

	`<SNL>  <PGA>  <PGV>  <PGD>  <PA3>  <PV3>  <PD3>  <TauC3>  <PGA Leading>  <PGV Leading>  <Epc. Dist.>  <S/N Ratio>`

- With station coordinate:

	`<SNL>  <PGA>  <PGV>  <PGD>  <PA3>  <PV3>  <PD3>  <TauC3>  <PGA Leading>  <PGV Leading>  <Epc. Dist.>  <S/N Ratio>  <Latitude>  <Longitude>  <Elevation>`

## License

Copyright&copy; 2019-2024 Benjamin Ming Yang

Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at

[Apache License 2.0](http://www.apache.org/licenses/LICENSE-2.0)

Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.