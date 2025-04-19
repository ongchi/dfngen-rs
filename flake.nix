{
  inputs = {
    fenix = {
      url = "github:nix-community/fenix";
      inputs.nixpkgs.follows = "nixpkgs";
    };
    flake-utils.url = "github:numtide/flake-utils";
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-24.11";
    rust-overlay.url = "github:oxalica/rust-overlay";
  };

  outputs = {
    self,
    fenix,
    flake-utils,
    nixpkgs,
    rust-overlay,
  }:
    flake-utils.lib.eachDefaultSystem (system: let
      toolchain = fenix.packages.${system}.stable.toolchain;
      overlays = [ (import rust-overlay) ];
      pkgs = import nixpkgs {
        inherit system overlays;
      };
    in {
      devShells.default = pkgs.mkShell {
        buildInputs = [
          toolchain
          pkgs.libiconv
          pkgs.cargo-audit
        ];
      };
    });
}
