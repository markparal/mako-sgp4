# mako-sgp4-wasm

WebAssembly bindings for [mako-sgp4](../README.md), so the propagator can run in a browser. Not published to crates.io.

## Build

```bash
rustup target add wasm32-unknown-unknown
cargo install wasm-pack
wasm-pack build --target web --profile wasm-release
```

Run from this directory. Copy `pkg/mako_sgp4_wasm.js` and `pkg/mako_sgp4_wasm_bg.wasm` to your website.

## Usage

```js
import init, { Satellite } from './mako_sgp4_wasm.js';
await init();

const sat = new Satellite(tleText);            // TLE or OMM KVN
const state = sat.propagate(60);               // [x, y, z, vx, vy, vz] 60 min after epoch, TEME [km, km/s]
const track = sat.trackGeodetic(0, 1440, 1);   // [lat, lon, alt, ...] for one day [deg, deg, km]
```

See `src/lib.rs` for the full API.
