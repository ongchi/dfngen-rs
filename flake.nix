{
  inputs = {
    fenix = {
      url = "github:nix-community/fenix";
      inputs.nixpkgs.follows = "nixpkgs";
    };
    flake-utils.url = "github:numtide/flake-utils";
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-24.11";
  };

  outputs = {
    self,
    fenix,
    flake-utils,
    nixpkgs,
  }:
    flake-utils.lib.eachDefaultSystem (system: let
      toolchain = with fenix.packages.${system};
        combine [
          stable.cargo
          stable.clippy
          stable.rust-src
          stable.rustc
          stable.rustfmt
          rust-analyzer
        ];
      pkgs = nixpkgs.legacyPackages.${system};
    in {
      devShells.default = pkgs.mkShell {
        buildInputs = [
          toolchain
          pkgs.libiconv
        ];
      };
    });
}
