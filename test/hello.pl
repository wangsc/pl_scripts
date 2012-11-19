#!/usr/bin/perl -w
use Inline C => <<'END_C';

void greet() {
  printf("Hello, world!\n");
}
END_C

greet;